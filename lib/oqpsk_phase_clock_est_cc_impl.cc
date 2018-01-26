/* -*- c++ -*- */
/* 
 * Copyright 2018 <+YOU OR YOUR COMPANY+>.
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#define _USE_MATH_DEFINES
#define NL0 1000

#include <cstdio>
#include <cmath>
#include <algorithm>
#include <functional>
#include <numeric>
#include <gnuradio/io_signature.h>
#include "oqpsk_phase_clock_est_cc_impl.h"
#include <math.h>
#include <volk/volk.h>


namespace gr {
  namespace oqpsk {

    oqpsk_phase_clock_est_cc::sptr
    oqpsk_phase_clock_est_cc::make()
    {
      return gnuradio::get_initial_sptr
        (new oqpsk_phase_clock_est_cc_impl());
    }

    /*
     * The private constructor
     */
    oqpsk_phase_clock_est_cc_impl::oqpsk_phase_clock_est_cc_impl()
      : gr::block("oqpsk_phase_clock_est_cc",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(1, 1, sizeof(gr_complex)))
    {
      gr_complex *top_b4_filt = (gr_complex *)volk_malloc(sizeof(gr_complex)*NL0, alignment);
      gr_complex *bot_b4_filt = (gr_complex *)volk_malloc(sizeof(gr_complex)*NL0, alignment);
      gr_complex *top_aft_filt = (gr_complex *)volk_malloc(sizeof(gr_complex)*NL0, alignment);
      gr_complex *bot_aft_filt = (gr_complex *)volk_malloc(sizeof(gr_complex)*NL0, alignment);
      gr_complex *top_b4_sum = (gr_complex *)volk_malloc(sizeof(gr_complex)*NL0, alignment);
      gr_complex *bot_b4_sum = (gr_complex *)volk_malloc(sizeof(gr_complex)*NL0, alignment);
      
      set_history(d_taps_q + d_sps);
      
      set_output_multiple(1);
    }

    /*
     * Our virtual destructor.
     */
    oqpsk_phase_clock_est_cc_impl::~oqpsk_phase_clock_est_cc_impl()
    {
      volk_free(top_b4_filt);
      volk_free(bot_b4_filt);
      volk_free(top_aft_filt);
      volk_free(bot_aft_filt);
      volk_free(top_b4_sum);
      volk_free(bot_b4_sum);
    }

    void
    oqpsk_phase_clock_est_cc_impl::forecast (int noutput_items, 
                                             gr_vector_int &ninput_items_required)
    {
      /* <+forecast+> e.g. ninput_items_required[0] = noutput_items */
    }

    int
    oqpsk_phase_clock_est_cc_impl::general_work (int noutput_items,
                                                 gr_vector_int &ninput_items,
                                                 gr_vector_const_void_star &input_items,
                                                 gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex*) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];
      float *out1 = (float *) output_items[1];
      int noi = ninput_items[0];

      // Generating input for filter q(t)
      gr_complex sample_phase_top, sample_phase_bot;
      float phase_arg;
      int nsamples = std::min(NL0, noi);
      for (int i = 0; i <= nsamples - 1; i++) {
        phase_arg = M_PI * i / noi;
        sample_phase_top = gr_expj(-phase_arg);
        sample_phase_bot = gr_expj(phase_arg);
        top_b4_filt[i] = in * sample_phase_top;
        bot_b4_filt[i] = in * sample_phase_bot;
      }
      
      d_q_filter->filter(nsamples, top_b4_filt, top_aft_filt);
      d_q_filter->filter(nsamples, bot_b4_filt, bot_aft_filt);
      
      gr_complex X, Y;
      
      X = std::inner_product(top_aft_filt, top_aft_filt + nsamples, top_b4_filt,
                             top_b4_sum, 0);
      Y = std::inner_product(bot_aft_filt, bot_aft_filt + nsamples, bot_b4_filt,
                             bot_b4_sum, 0);
      
      d_phase = 1 / 4 * (atan2(X.imag(), X.real()) + atan2(Y.imag(), Y.real()));
      d_delay = 1 / (4 * M_PI) * (-atan2(X.imag(), X.real()) + atan2(Y.imag(), Y.real()));
      
      gr_complex est_phase;
      est_phase = gr_expj(-d_phase);
      std::transform(in, in + nsamples, out, 
                     std::bind1st(std::multiplies<gr_complex>(),est_phase));
      *out1 = d_phase;
      *out2 = d_delay;
      produce(0, nsamples);
      produce(1, 1);
      produce(2, 1);
      
      
      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each (nsamples);

      // Tell runtime system how many output items we produced.
      return WORK_CALLED_PRODUCE;
    }

  } /* namespace oqpsk */
} /* namespace gr */

