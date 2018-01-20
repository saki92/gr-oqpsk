/* -*- c++ -*- */

#define OQPSK_API

%include "gnuradio.i"			// the common stuff

//load generated python docstrings
%include "oqpsk_swig_doc.i"

%{
#include "oqpsk/pfb_oqpsk_clock_sync_cc.h"
#include "oqpsk/oqpsk_phase_sync_cc.h"
%}


%include "oqpsk/pfb_oqpsk_clock_sync_cc.h"
GR_SWIG_BLOCK_MAGIC2(oqpsk, pfb_oqpsk_clock_sync_cc);
%include "oqpsk/oqpsk_phase_sync_cc.h"
GR_SWIG_BLOCK_MAGIC2(oqpsk, oqpsk_phase_sync_cc);
