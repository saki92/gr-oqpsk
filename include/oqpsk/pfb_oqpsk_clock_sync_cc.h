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


#ifndef INCLUDED_OQPSK_PFB_OQPSK_CLOCK_SYNC_CC_H
#define INCLUDED_OQPSK_PFB_OQPSK_CLOCK_SYNC_CC_H

#include <oqpsk/api.h>
#include <gnuradio/filter/fir_filter.h>
#include <gnuradio/block.h>

namespace gr {
  namespace oqpsk {

    /*!
     * \brief <+description of block+>
     * \ingroup oqpsk
     *
     */
    class OQPSK_API pfb_oqpsk_clock_sync_cc : virtual public gr::block
    {
     public:
      typedef boost::shared_ptr<pfb_oqpsk_clock_sync_cc> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of oqpsk::pfb_oqpsk_clock_sync_cc.
       *
       * To avoid accidental use of raw pointers, oqpsk::pfb_oqpsk_clock_sync_cc's
       * constructor is in a private implementation
       * class. oqpsk::pfb_oqpsk_clock_sync_cc::make is the public interface for
       * creating new instances.
       */
      static sptr make(double sps, float loop_bw, 
                       const std::vector<float> &taps, 
                       unsigned int filter_size, 
                       float init_phase, float max_rate_deviation);
                       
      virtual void update_gains() = 0;

      /*!
       * Resets the filterbank's filter taps with the new prototype filter.
       */
      virtual void update_taps(const std::vector<float> &taps) = 0;

      /*!
       * Used to set the taps of the filters in the filterbank and
       * differential filterbank.
       *
       * WARNING: this should not be used externally and will be moved
       * to a private function in the next API.
       */
      virtual void set_taps(const std::vector<float> &taps,
			    std::vector< std::vector<float> > &ourtaps,
			    std::vector<gr::filter::kernel::fir_filter_fff*> &ourfilter) = 0;

      /*!
       * Returns all of the taps of the matched filter
       */
      virtual std::vector< std::vector<float> > taps() const = 0;

      /*!
       * Returns all of the taps of the derivative filter
       */
      virtual std::vector< std::vector<float> > diff_taps() const = 0;

      /*!
       * Returns the taps of the matched filter for a particular channel
       */
      virtual std::vector<float> channel_taps(int channel) const = 0;

      /*!
       * Returns the taps in the derivative filter for a particular channel
       */
      virtual std::vector<float> diff_channel_taps(int channel) const = 0;

      /*!
       * Return the taps as a formatted string for printing
       */
      virtual std::string taps_as_string() const = 0;

      /*!
       * Return the derivative filter taps as a formatted string for printing
       */
      virtual std::string diff_taps_as_string() const = 0;


      /*******************************************************************
       SET FUNCTIONS
      *******************************************************************/

      /*!
       * \brief Set the loop bandwidth
       *
       * Set the loop filter's bandwidth to \p bw. This should be
       * between 2*pi/200 and 2*pi/100 (in rads/samp). It must also be
       * a positive number.
       *
       * When a new damping factor is set, the gains, alpha and beta,
       * of the loop are recalculated by a call to update_gains().
       *
       * \param bw    (float) new bandwidth
       */
      virtual void set_loop_bandwidth(float bw) = 0;

      /*!
       * \brief Set the loop damping factor
       *
       * Set the loop filter's damping factor to \p df. The damping
       * factor should be sqrt(2)/2.0 for critically damped systems.
       * Set it to anything else only if you know what you are
       * doing. It must be a number between 0 and 1.
       *
       * When a new damping factor is set, the gains, alpha and beta,
       * of the loop are recalculated by a call to update_gains().
       *
       * \param df    (float) new damping factor
       */
      virtual void set_damping_factor(float df) = 0;

      /*!
       * \brief Set the loop gain alpha
       *
       * Set's the loop filter's alpha gain parameter.
       *
       * This value should really only be set by adjusting the loop
       * bandwidth and damping factor.
       *
       * \param alpha    (float) new alpha gain
       */
      virtual void set_alpha(float alpha) = 0;

      /*!
       * \brief Set the loop gain beta
       *
       * Set's the loop filter's beta gain parameter.
       *
       * This value should really only be set by adjusting the loop
       * bandwidth and damping factor.
       *
       * \param beta    (float) new beta gain
       */
      virtual void set_beta(float beta) = 0;

      /*!
       * Set the maximum deviation from 0 d_rate can have
       */
      virtual void set_max_rate_deviation(float m) = 0;

      /*******************************************************************
       GET FUNCTIONS
      *******************************************************************/

      /*!
       * \brief Returns the loop bandwidth
       */
      virtual float loop_bandwidth() const = 0;

      /*!
       * \brief Returns the loop damping factor
       */
      virtual float damping_factor() const = 0;

      /*!
       * \brief Returns the loop gain alpha
       */
      virtual float alpha() const = 0;

      /*!
       * \brief Returns the loop gain beta
       */
      virtual float beta() const = 0;

      /*!
       * \brief Returns the current clock rate
       */
      virtual float clock_rate() const = 0;

      /*!
       * \brief Returns the current error of the control loop.
       */
      virtual float error() const = 0;

      /*!
       * \brief Returns the current rate of the control loop.
       */
      virtual float rate() const = 0;

      /*!
       * \brief Returns the current phase arm of the control loop.
       */
      virtual float phase() const = 0;
    };

  } // namespace oqpsk
} // namespace gr

#endif /* INCLUDED_OQPSK_PFB_OQPSK_CLOCK_SYNC_CC_H */

