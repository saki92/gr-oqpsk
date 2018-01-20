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

#include <gnuradio/io_signature.h>
#include "oqpsk_phase_sync_cc_impl.h"
#include <gnuradio/expj.h>
#include <gnuradio/sincos.h>
#include <gnuradio/math.h>

namespace gr {
  namespace oqpsk {

    oqpsk_phase_sync_cc::sptr
    oqpsk_phase_sync_cc::make(float loop_bw)
    {
      return gnuradio::get_initial_sptr
        (new oqpsk_phase_sync_cc_impl(loop_bw));
    }

    static int ios[] = { sizeof(gr_complex), sizeof(float), sizeof(float), sizeof(float) };
    static std::vector<int> iosig(ios, ios+sizeof(ios)/sizeof(int));
    oqpsk_phase_sync_cc_impl::oqpsk_phase_sync_cc_impl(float loop_bw)
      : gr::sync_block("oqpsk_phase_sync_cc",
              gr::io_signature::make(1, 2, sizeof(gr_complex)),
              gr::io_signature::makev(1, 4, iosig)),
        blocks::control_loop(loop_bw, 1.0, -1.0),
	d_error(0), d_noise(1.0)
    {}

    oqpsk_phase_sync_cc_impl::~oqpsk_phase_sync_cc_impl()
    {
    }

    float
    oqpsk_phase_sync_cc_impl::phase_estimator(gr_complex sample,
                                     gr_complex co_sample)
    {
      float K = 2 * abs(sample) / d_noise;
      return ((blocks::tanhf_lut(K * sample.imag()) * co_sample.real()) -
              (blocks::tanhf_lut(K * sample.real()) * co_sample.imag()));
    } //tanh from volk can be used to further reduce quantization noise

    int
    oqpsk_phase_sync_cc_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      const gr_complex *iptr = (gr_complex *) input_items[0];
      const gr_complex *iptr1 = (gr_complex *) input_items[1];
      gr_complex *optr = (gr_complex *) output_items[0];
      float *freq_optr  = output_items.size() >= 2 ? (float *) output_items[1] : NULL;
      float *phase_optr = output_items.size() >= 3 ? (float *) output_items[2] : NULL;
      float *error_optr = output_items.size() >= 4 ? (float *) output_items[3] : NULL;
      
      gr_complex nco_out;

      for(int i = 0; i < noutput_items; i++) {
        nco_out = gr_expj(-d_phase);
        optr[i] = iptr[i] * nco_out;

        d_error = phase_estimator(optr[i], iptr1[i] * nco_out);
        d_error = gr::branchless_clip(d_error, 1.0);

        advance_loop(d_error);
        phase_wrap();
        frequency_limit();

        if (freq_optr != NULL)
          freq_optr[i] = d_freq;
        if (phase_optr != NULL)
          phase_optr[i] = d_phase;
        if (error_optr != NULL)
          error_optr[i] = d_error;
      }
      return noutput_items;
    }

  } /* namespace oqpsk */
} /* namespace gr */

