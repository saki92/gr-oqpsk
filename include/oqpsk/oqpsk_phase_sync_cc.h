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


#ifndef INCLUDED_OQPSK_OQPSK_PHASE_SYNC_CC_H
#define INCLUDED_OQPSK_OQPSK_PHASE_SYNC_CC_H

#include <oqpsk/api.h>
#include <gnuradio/sync_block.h>
#include <gnuradio/blocks/control_loop.h>

namespace gr {
  namespace oqpsk {

    /*!
     * \brief <+description of block+>
     * \ingroup oqpsk
     *
     */
    class OQPSK_API oqpsk_phase_sync_cc : virtual public sync_block, virtual public blocks::control_loop
    {
     public:
      typedef boost::shared_ptr<oqpsk_phase_sync_cc> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of oqpsk::oqpsk_phase_sync_cc.
       *
       * To avoid accidental use of raw pointers, oqpsk::oqpsk_phase_sync_cc's
       * constructor is in a private implementation
       * class. oqpsk::oqpsk_phase_sync_cc::make is the public interface for
       * creating new instances.
       */
      static sptr make(float loop_bw);
    };

  } // namespace oqpsk
} // namespace gr

#endif /* INCLUDED_OQPSK_OQPSK_PHASE_SYNC_CC_H */

