/*
 *  gpe.cpp
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

 #include "gpe.h"

 // C++ includes:
 #include <limits>

 // Includes from libnestutil:
 #include "numerics.h"

 // Includes from nestkernel:
 #include "event_delivery_manager_impl.h"
 #include "exceptions.h"
 #include "kernel_manager.h"
 #include "universal_data_logger_impl.h"

 // Includes from sli:
 #include "dict.h"
 #include "dictutils.h"
 #include "doubledatum.h"
 #include "integerdatum.h"

 /* ----------------------------------------------------------------
  * Recordables map
  * ---------------------------------------------------------------- */

 nest::RecordablesMap< nest::GPe > nest::GPe::recordablesMap_;

 namespace nest
 {
 // Override the create() method with one call to RecordablesMap::insert_()
 // for each quantity to be recorded.
 template <>
 void
 RecordablesMap< GPe >::create()
 {
   // use standard names whereever you can for consistency!
   insert_( names::V_m, &GPe::get_V_m_ );
   insert_( names::U_m, &GPe::get_U_m_ );
 }
 }


 /* ----------------------------------------------------------------
  * Default constructors defining default parameters and state
  * ---------------------------------------------------------------- */

 nest::GPe::Parameters_::Parameters_()
   : a_( 0.03 )                                      // a
   , b_( -2.0 )                                      // b
   , c_( -50.0 )                                     // c without unit
   , d_( 100.0 )                                     // d
   , I_e_( 0.0 )                                     // pA
   , V_th_( 35.0 )                                   // mV
   , V_min_( -std::numeric_limits< double >::max() ) // mV
   , Kplus_( 0.7 )                                   // k without unit
   , V_m_( -60.0 )                                   // mV
   , V_T_( -40.0 )                                   // mV
   , C_m_( 100.0 )                                   // pF

   , c_1_( 0.5 )                                  // no unit
   , c_2_( 0.5 )                                   // no unit
   , y1_( 0.3 )                                      // no unit
   , y2_( 0.3 )                                      // no unit

   , tau_syn_( 1, 2.0 ) // ms
   , E_rev_( 1, 0.0 )
 {
 }

 nest::GPe::State_::State_()
   : v_( -65.0 ) // membrane potential
   , u_( 0.0 )   // membrane recovery variable
   , I_( 0.0 )   // input current

   , y_(1, 0.0)
 {
 }

 /* ----------------------------------------------------------------
  * Parameter and state extractions and manipulation functions
  * ---------------------------------------------------------------- */

 void
 nest::GPe::Parameters_::get( DictionaryDatum& d ) const
 {
   def< double >( d, names::I_e, I_e_ );
   def< double >( d, names::V_th, V_th_ ); // threshold value
   def< double >( d, names::V_min, V_min_ );
   def< double >( d, names::a, a_ );
   def< double >( d, names::b, b_ );
   def< double >( d, names::c, c_ );
   def< double >( d, names::d, d_ );
   def< double >( d, names::Kplus, Kplus_ );
   def< double >( d, names::V_m, V_m_ );
   def< double >( d, names::V_T, V_T_ );
   def< double >( d, names::C_m, C_m_ );
   def< double >( d, names::c_1, c_1_ );
   def< double >( d, names::c_2, c_2_ );
   def< double >( d, names::y1, y1_ );
   def< double >( d, names::y2, y2_ );

   def< size_t >( d, names::n_receptors, n_receptors() );
   ArrayDatum E_rev_ad( E_rev_ );
   def< ArrayDatum >( d, names::E_rev, E_rev_ad );
   ArrayDatum tau_syn_ad( tau_syn_ );
   def< ArrayDatum >( d, names::tau_syn, tau_syn_ad );
 }

 void
 nest::GPe::Parameters_::set( const DictionaryDatum& d )
 {

   updateValue< double >( d, names::V_th, V_th_ );
   updateValue< double >( d, names::V_min, V_min_ );
   updateValue< double >( d, names::I_e, I_e_ );
   updateValue< double >( d, names::a, a_ );
   updateValue< double >( d, names::b, b_ );
   updateValue< double >( d, names::c, c_ );
   updateValue< double >( d, names::d, d_ );

   updateValue< double >( d, names::Kplus, Kplus_ );
   updateValue< double >( d, names::V_m, V_m_ );
   updateValue< double >( d, names::V_T, V_T_ );
   updateValue< double >( d, names::C_m, C_m_ );

   updateValue< double >( d, names::c_1, c_1_ );
   updateValue< double >( d, names::c_2, c_2_ );
   updateValue< double >( d, names::y1, y1_ );
   updateValue< double >( d, names::y2, y2_ );

   const size_t old_n_receptors = n_receptors();
   bool Erev_flag =
     updateValue< std::vector< double > >( d, names::E_rev, E_rev_ );
   bool tau_flag =
     updateValue< std::vector< double > >( d, names::tau_syn, tau_syn_ );

   // receptor arrays have been modified
   if ( Erev_flag || tau_flag )
   {
     if ( ( E_rev_.size() != old_n_receptors
            || tau_syn_.size() != old_n_receptors )
       and ( not Erev_flag || not tau_flag ) )
     {
       throw BadProperty(
         "If the number of receptor ports is changed, both arrays "
         "E_rev and tau_syn must be provided." );
     }
     if ( E_rev_.size() != tau_syn_.size() )
     {
       throw BadProperty(
         "The reversal potential, and synaptic time constant arrays "
         "must have the same size." );
     }
     if ( tau_syn_.size() < old_n_receptors && has_connections_ )
     {
       throw BadProperty(
         "The neuron has connections, therefore the number of ports cannot be "
         "reduced." );
     }
     for ( size_t i = 0; i < tau_syn_.size(); ++i )
     {
       if ( tau_syn_[ i ] <= 0 )
       {
         throw BadProperty(
           "All synaptic time constants must be strictly positive" );
       }
     }
   }
 }

 void
 nest::GPe::State_::get( DictionaryDatum& d, const Parameters_& ) const
 {
   def< double >( d, names::U_m, u_ ); // Membrane potential recovery variable
   def< double >( d, names::V_m, v_ ); // Membrane potential

   std::vector< double >* g = new std::vector< double >();
   for (size_t i=0; i<y_.size(); ++i)
     g->push_back(y_[i]);
   ( *d )[ names::g ] = DoubleVectorDatum( g );
 }

 void
 nest::GPe::State_::set( const DictionaryDatum& d, const Parameters_& )
 {
   updateValue< double >( d, names::U_m, u_ );
   updateValue< double >( d, names::V_m, v_ );
 }

 nest::GPe::Buffers_::Buffers_( GPe& n )
   : logger_( n )
 {
 }

 nest::GPe::Buffers_::Buffers_( const Buffers_&, GPe& n )
   : logger_( n )
 {
 }


 /* ----------------------------------------------------------------
  * Default and copy constructor for node
  * ---------------------------------------------------------------- */

 nest::GPe::GPe()
   : Archiving_Node()
   , P_()
   , S_()
   , B_( *this )
 {
   recordablesMap_.create();
 }

 nest::GPe::GPe( const GPe& n )
   : Archiving_Node( n )
   , P_( n.P_ )
   , S_( n.S_ )
   , B_( n.B_, *this )
 {
 }


 /* ----------------------------------------------------------------
  * Node initialization functions
  * ---------------------------------------------------------------- */

 void
 nest::GPe::init_state_( const Node& proto )
 {
   const GPe& pr = downcast< GPe >( proto );
   S_ = pr.S_;
 }

 void
 nest::GPe::init_buffers_()
 {
   B_.spikes_.clear();   // includes resize
   B_.currents_.clear(); // includes resize
   B_.logger_.reset();   // includes resize
   Archiving_Node::clear_history();
 }

 void
 nest::GPe::calibrate()
 {
   B_.spikes_.resize( P_.n_receptors() );
   S_.y_.resize(P_.n_receptors(), 0.0);

   B_.logger_.init();
 }


 /* ----------------------------------------------------------------
  * Update and spike handling functions
  */

 void
 nest::GPe::update( Time const& origin, const long from, const long to )
 {
   assert(
     to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
   assert( from < to );

   const double h = Time::get_resolution().get_ms();

   for ( long lag = from; lag < to; ++lag )
   {
     // neuron is never refractory

     // Get the value current of the synapses
     // This models adds this state variables (S_):
     // - y_[i] is the contribution to membrane potential of receptor i.
     // and this parameters (P_):
     // - tau_syn_[i] is the temporal constant of receptor i.
     // - E_rev_[i] is the receptor reversal potential


/*
     double I_syn = 0.0;
     for (size_t i=0; i<P_.n_receptors(); i++) {
       S_.y_[i] -= S_.y_[i] / P_.tau_syn_[i];
       S_.y_[i] += B_.spikes_[i].get_value( lag );
       I_syn += S_.y_[i] * (P_.E_rev_[i] - S_.v_);
     }
*/
     // Unrolling loop for dealing with receptors modulated by dopamine
     double I_syn = 0.0;
     // AMPA (modulated by dopamine)
     S_.y_[0] -= h * (S_.y_[0] / P_.tau_syn_[0]);
     S_.y_[0] += B_.spikes_[0].get_value( lag );
     I_syn += S_.y_[0] * (P_.E_rev_[0] - S_.v_) * (1 - P_.c_1_ * P_.y2_);
     // NMDA (modulated by dopamine)
     S_.y_[1] -= h * (S_.y_[1] / P_.tau_syn_[1]);
     S_.y_[1] += B_.spikes_[1].get_value( lag );
//     double B = 1.0/(1.0+(0.28)*exp(-0.062*S_.v_));
//     I_syn += S_.y_[1] * (P_.E_rev_[1] - S_.v_) * (1 - P_.c_1_ * P_.y2_) * B;
     I_syn += S_.y_[1] * (P_.E_rev_[1] - S_.v_) * (1 - P_.c_1_ * P_.y2_);
     // GABAa (modulated by dopamine)
     S_.y_[2] -= h * (S_.y_[2] / P_.tau_syn_[2]);
     S_.y_[2] += B_.spikes_[2].get_value( lag );
     I_syn += S_.y_[2] * (P_.E_rev_[2] - S_.v_) * (1 - P_.c_2_ * P_.y2_);
     // GABAa (modulated by dopamine)
     S_.y_[3] -= h * (S_.y_[3] / P_.tau_syn_[3]);
     S_.y_[3] += B_.spikes_[3].get_value( lag );
     I_syn += S_.y_[3] * (P_.E_rev_[3] - S_.v_) * (1 - P_.c_2_ * P_.y2_);


     // use numerics published in Izhikevich (2007)

     S_.v_ += h * (P_.Kplus_ * (S_.v_ - P_.V_m_) * (S_.v_ - P_.V_T_) - S_.u_ + S_.I_
                   + P_.I_e_ + I_syn ) / P_.C_m_;
     S_.u_ += h * P_.a_ * ( P_.b_ * (S_.v_ - P_.V_m_) - S_.u_ );


     // lower bound of membrane potential
     S_.v_ = ( S_.v_ < P_.V_min_ ? P_.V_min_ : S_.v_ );

     // threshold crossing
     if ( S_.v_ >= P_.V_th_ )
     {
       S_.v_ = P_.c_;
       S_.u_ = S_.u_ + P_.d_;

       // compute spike time
       set_spiketime( Time::step( origin.get_steps() + lag + 1 ) );

       SpikeEvent se;
       kernel().event_delivery_manager.send( *this, se, lag );
     }

     // set new input current
     S_.I_ = B_.currents_.get_value( lag );

     // voltage logging
     B_.logger_.record_data( origin.get_steps() + lag );
   }
 }



void
nest::GPe::handle( SpikeEvent& e )
{
 assert( e.get_delay_steps() > 0 );
 assert( ( e.get_rport() > 0 ) && ( ( size_t ) e.get_rport() <= P_.n_receptors() ) );

 B_.spikes_[ e.get_rport() - 1 ].add_value(
   e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ),
   e.get_weight() * e.get_multiplicity() );
}

void
nest::GPe::handle( CurrentEvent& e )
{
 assert( e.get_delay_steps() > 0 );

 const double c = e.get_current();
 const double w = e.get_weight();
 B_.currents_.add_value(
   e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ),
   w * c );
}

void
nest::GPe::handle( DataLoggingRequest& e )
{
 B_.logger_.handle( e );
}
