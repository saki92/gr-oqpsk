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

#ifndef INCLUDED_OQPSK_PFB_OQPSK_CLOCK_SYNC_CC_IMPL_H
#define INCLUDED_OQPSK_PFB_OQPSK_CLOCK_SYNC_CC_IMPL_H

#include <oqpsk/pfb_oqpsk_clock_sync_cc.h>

using namespace gr::filter;

namespace gr {
  namespace oqpsk {

    class pfb_oqpsk_clock_sync_cc_impl : public pfb_oqpsk_clock_sync_cc
    {
     private:
      bool   d_updated;
      double d_sps;
      double d_sample_num;
      float  d_loop_bw;
      float  d_damping;
      float  d_alpha;
      float  d_beta;

      int                                  d_nfilters;
      int                                  d_taps_per_filter;
      std::vector<kernel::fir_filter_fff*> d_filters;
      std::vector<kernel::fir_filter_fff*> d_diff_filters;
      std::vector< std::vector<float> >    d_taps;
      std::vector< std::vector<float> >    d_dtaps;
      std::vector<float>                   d_updated_taps;
      kernel::fir_filter_fff              *d_diff_filter_one;

      float d_k_r;
      float d_k_i;
      float d_rate;
      float d_rate_r;
      float d_rate_i;
      float d_rate_r_i;
      float d_rate_r_f;
      float d_rate_i_i;
      float d_rate_i_f;
      float d_max_dev;
      int   d_filtnum_r;
      int   d_filtnum_i;
      int   d_osps;
      float d_error_r;
      float d_error_i;
      int   d_out_idx;
      int   d_count_dist;

      uint64_t d_old_in, d_new_in, d_last_out;

      void create_diff_taps(const std::vector<float> &newtaps,
			    std::vector<float> &difftaps);

     public:
      pfb_oqpsk_clock_sync_cc_impl(double sps, float loop_bw, 
      			const std::vector<float> &taps, 
      			unsigned int filter_size = 32, 
      			float init_phase = 0, 
      			float max_rate_deviation = 1.5);
      			
      ~pfb_oqpsk_clock_sync_cc_impl();

      void setup_rpc();

      void update_gains();

      void forecast(int noutput_items, gr_vector_int &ninput_items_required);

      void update_taps(const std::vector<float> &taps);

      void set_taps(const std::vector<float> &taps,
		    std::vector< std::vector<float> > &ourtaps,
		    std::vector<kernel::fir_filter_fff*> &ourfilter);

      std::vector< std::vector<float> > taps() const;
      std::vector< std::vector<float> > diff_taps() const;
      std::vector<float> channel_taps(int channel) const;
      std::vector<float> diff_channel_taps(int channel) const;
      std::string taps_as_string() const;
      std::string diff_taps_as_string() const;

      void set_loop_bandwidth(float bw);
      void set_damping_factor(float df);
      void set_alpha(float alpha);
      void set_beta(float beta);
      void set_max_rate_deviation(float m)
      {
	d_max_dev = m;
      }

      float loop_bandwidth() const;
      float damping_factor() const;
      float alpha() const;
      float beta() const;
      float clock_rate() const;

      float error() const;
      float rate() const;
      float phase() const;

      /*******************************************************************
       *******************************************************************/

      bool check_topology(int ninputs, int noutputs);

      int general_work(int noutput_items,
		       gr_vector_int &ninput_items,
		       gr_vector_const_void_star &input_items,
		       gr_vector_void_star &output_items);
      /*
      int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);
      */
    };

  } // namespace oqpsk
} // namespace gr

#endif /* INCLUDED_OQPSK_PFB_OQPSK_CLOCK_SYNC_CC_IMPL_H */

