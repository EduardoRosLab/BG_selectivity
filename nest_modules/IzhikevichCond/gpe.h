/*
 *  gpe.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef GPe_H
#define GPe_H

// Includes from nestkernel:
#include "archiving_node.h"
#include "connection.h"
#include "event.h"
#include "nest_types.h"
#include "ring_buffer.h"
#include "universal_data_logger.h"

namespace nest
{
/* BeginDocumentation
   Name: izhikevich_complex - Izhikevich complex neuron model

   Description:
   Implementation of the spiking neuron model introduced by Izhikevich [1]. The
   dynamics are given by:

        C(dv/dt) = k * (v - v_r) * (v - v_t) - u + I
        du/dt = a * (b * (v - v_r) - u)

       if v >= V_th:
         v is set to c
         u is incremented by d

       v jumps on each spike arrival by the weight of the spike.

   As published in [1], the numerics differs from the standard forward Euler
   technique in two ways:
   1) the new value of u is calculated based on the new value of v, rather than
   the previous value
   2) the variable v is updated using a time step half the size of that used to
   update variable u.

   Parameters:
   The following parameters can be set in the status dictionary.

   V_m        double - Membrane potential in mV
   U_m        double - Membrane potential recovery variable
   V_th       double - Spike threshold in mV.
   I_e        double - Constant input current in pA. (R=1)
   V_min      double - Absolute lower value for the membrane potential.
   a          double - describes time scale of recovery variable
   b          double - sensitivity of recovery variable
   c          double - after-spike reset value of V_m
   d          double - after-spike reset value of U_m

   Kplus      double - Abstract parameter (k in [1])
   V_m        double - Resting membrane potential (v_r in [1])
   V_T        double - Instantaneous threshold potential (v_t in [1])
   C_m        double - Membrane capacitance (C in [1])

   Dopamine-related parameters
   c_1        double - beta_1 in [2]
   c_2        double - beta_2 in [2]
   y1         double - d1 in [2] represents the level of dopamine in the system
              affecting to D1 receptors.
   y2         double - d2 in [2] represents the level of dopamine in the system
              affecting to D2 receptors.




   References:
   [1] Izhikevich, E. M. (2007). Dynamical Systems in Neuroscience Computational
   Neuroscience. Dynamical Systems (Vol. 25).
   https://doi.org/10.1017/S0143385704000173
   [2] Fountas, Z. (2016). Action selection in the rhythmic brain: The role of
   the basal ganglia and tremor.

   Sends: SpikeEvent

   Receives: SpikeEvent, CurrentEvent, DataLoggingRequest
   FirstVersion: 2017
   Author: Álvaro González Redondo. Based on Hanuschkin, Morrison, Kunkel work.
   SeeAlso: izhikevich, iaf_psc_delta, mat2_psc_exp
*/

class GPe : public Archiving_Node
{

public:
  GPe();
  GPe( const GPe& );

  /**
   * Import sets of overloaded virtual functions.
   * @see Technical Issues / Virtual Functions: Overriding, Overloading, and
   * Hiding
   */
  using Node::handle;
  using Node::handles_test_event;

  void handle( DataLoggingRequest& );
  void handle( SpikeEvent& );
  void handle( CurrentEvent& );

  port handles_test_event( DataLoggingRequest&, rport );
  port handles_test_event( SpikeEvent&, rport );
  port handles_test_event( CurrentEvent&, rport );

  port send_test_event( Node&, rport, synindex, bool );

  void get_status( DictionaryDatum& ) const;
  void set_status( const DictionaryDatum& );

private:
  friend class RecordablesMap< GPe >;
  friend class UniversalDataLogger< GPe >;

  void init_state_( const Node& proto );
  void init_buffers_();
  void calibrate();

  void update( Time const&, const long, const long );

  // ----------------------------------------------------------------

  /**
   * Independent parameters of the model.
   */
  struct Parameters_
  {
    double a_;
    double b_;
    double c_;
    double d_;

    /** External DC current */
    double I_e_;

    /** Threshold */
    double V_th_;

    /** Lower bound */
    double V_min_;

    /** Abstract parameter **/
    double Kplus_;

    /** Resting membrane potential **/
    double V_m_;

    /** Instantaneous threshold potential **/
    double V_T_;

    /** Membrane capacitance **/
    double C_m_;

    /** Betas (parameters of dopamine modulation) **/
    double c_1_, c_2_, c_3_;

    /** Dopamine levels **/
    double y1_, y2_;

    /** Time constants of synaptic currents in ms */
    std::vector< double > tau_syn_;

    /** reversal potentials in mV */
    std::vector< double > E_rev_;

    /** boolean flag which indicates whether the neuron has connections */
    bool has_connections_;

    Parameters_(); //!< Sets default parameter values

    void get( DictionaryDatum& ) const; //!< Store current values in dictionary
    void set( const DictionaryDatum& ); //!< Set values from dicitonary

    //! Return the number of receptor ports
    inline size_t n_receptors() const { return E_rev_.size(); }
  };

  // ----------------------------------------------------------------

  /**
   * State variables of the model.
   */
  struct State_
  {
    double v_; // membrane potential
    double u_; // membrane recovery variable
    double I_; // input current

    std::vector< double > y_;        //!< neuron state

    State_(); //!< Default initialization

    void get( DictionaryDatum&, const Parameters_& ) const;
    void set( const DictionaryDatum&, const Parameters_& );
  };

  // ----------------------------------------------------------------

  /**
   * Buffers of the model.
   */
  struct Buffers_
  {
    /**
     * Buffer for recording
     */
    Buffers_( GPe& );
    Buffers_( const Buffers_&, GPe& );
    UniversalDataLogger< GPe > logger_;

    /** buffers and sums up incoming spikes/currents */
    std::vector< RingBuffer > spikes_;
    RingBuffer currents_;
  };

  // ----------------------------------------------------------------

  /**
   * Internal variables of the model.
   */
  struct Variables_
  {
  };

  // Access functions for UniversalDataLogger -----------------------

  //! Read out the membrane potential
  double
  get_V_m_() const
  {
    return S_.v_;
  }
  //! Read out the recovery variable
  double
  get_U_m_() const
  {
    return S_.u_;
  }

  // ----------------------------------------------------------------

  Parameters_ P_;
  State_ S_;
  Variables_ V_;
  Buffers_ B_;

  //! Mapping of recordables names to access functions
  static RecordablesMap< GPe > recordablesMap_;
  /** @} */
};


inline port
GPe::send_test_event( Node& target, rport receptor_type, synindex, bool )
{
  SpikeEvent e;
  e.set_sender( *this );

  return target.handles_test_event( e, receptor_type );
}

inline port
GPe::handles_test_event( SpikeEvent&, rport receptor_type )
{
  if ( receptor_type <= 0
    || receptor_type > static_cast< port >( P_.n_receptors() ) )
  {
    throw IncompatibleReceptorType( receptor_type, get_name(), "SpikeEvent" );
  }

  P_.has_connections_ = true;
  return receptor_type;
}

inline port
GPe::handles_test_event( CurrentEvent&, rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}

inline port
GPe::handles_test_event( DataLoggingRequest& dlr, rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return B_.logger_.connect_logging_device( dlr, recordablesMap_ );
}

inline void
GPe::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  S_.get( d, P_ );
  Archiving_Node::get_status( d );
  ( *d )[ names::recordables ] = recordablesMap_.get_list();
}

inline void
GPe::set_status( const DictionaryDatum& d )
{
  Parameters_ ptmp = P_; // temporary copy in case of errors
  ptmp.set( d );         // throws if BadProperty
  State_ stmp = S_;      // temporary copy in case of errors
  stmp.set( d, ptmp );   // throws if BadProperty

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  Archiving_Node::set_status( d );

  // if we get here, temporaries contain consistent set of properties
  P_ = ptmp;
  S_ = stmp;
}

} // namespace nest

#endif /* #ifndef GPe_H */
