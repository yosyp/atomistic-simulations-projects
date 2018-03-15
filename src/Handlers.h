/*!
* @file Handlers.h
* @author Maxim Shugaev
* @date 22 Jan 2018
* @brief This file contains a set of classes used for calculation of system properties and controlling the simulation.
*  It includes calculation of Energy, Temperature, and Pressure as well as assignment of initial velocity and using Berendsen thermostat and barostat. */

#pragma once
#include "System.h"
#include <vector>

namespace MSE6270_MD {
	/*! An interface that declares methods used for calculation and controlling simulation properties */
	class Handler_Base {
	public:
		Handler_Base();												//!< Initialize Handler_Base
		virtual ~Handler_Base() = default;							//destructor must be virtual if inheritance is used
		void update(System &system);								//!< Update system
		void activate();											//!< Activate handler
		void deactivate();											//!< Deactivate handler
		void set_update_frequency(int n);							//!< Set the frequency (every "n"'s step) how often the system is updated
		virtual std::string info() const;							//!< Return a string with information about a handler				
	private:
		virtual void update_impl(System &system) = 0;				//!< This function implements a method updating the system. Override this function to implement a particular algorithm

		bool active;												/*!< A flag that determines if the handler is active */											
		int update_n;												/*!< The frequency of system update */
	};

	/*! This class is responsible for assignment of the velocity to particles according to the Maxwell distribution. 3KT energy is assigned per particle, 
		and half of it will be transferred to the potential energy (I assume that the initial configuration is quenched). More about initialization of velocities in Frenkel & Smit pp 56-57 */
	class Handler_Velocity : public Handler_Base {					//this class implements Handler_Base interface
	public:
		Handler_Velocity(double T0 = 300.0);						//Create a handler assigning T0 by setting particle velocities at the 0-th time step
		virtual ~Handler_Velocity() = default;
		virtual std::string info() const override;					//!< Return a string with information about a handler
	protected:
		virtual void update_impl(System &system) override;			//!< Assign velocities at the 0-th step

		double T0;													/*!< The temperature assigned by setting particle velocities */
	};

	/*! This class is responsible for calculation of the total potential and kinetics energies of the system */
	class Handler_Energy : public Handler_Base {					//this class implements Handler_Base interface
	public:
		Handler_Energy();
		virtual ~Handler_Energy() = default;
		double get_Ep() const;										//!< Return the total potential energy of the system
		double get_Ek() const;										//!< Return the total kinetic energy of the system
		double get_Etot() const;									//!< Return the total energy of the system
	protected:
		virtual void update_impl(System &system) override;			//!< Calculate Ep, Ek, and Etot

		double Ep;													/*!< The total potential energy of the system */
		double Ek;													/*!< The total kinetic energy of the system */
		double Etot;												/*!< The total energy of the system */
	};

	/*! This class is responsible for Temperature calculation */
	class Handler_Temperature : public Handler_Base {				//this class implements Handler_Base interface
	public:
		Handler_Temperature();
		virtual ~Handler_Temperature() = default;
		double get_T() const;										//!< Return the Temperature of the system
		void exclude_groups(const std::vector<int> &list);			//!< Set a list of groups (khist) that are excluded during temperature calculation
	protected:
		virtual void update_impl(System &system) override;			//!< Calculate the Temperature [K], estimated based on the average atomic kinetic energy: <Ek> = 1.5*kb*T
		bool exclude(int group) const;								//!< Check if a particular group in the exclusion list

		double T;													/*!< The Temperature of the system */
		std::vector<int> exclude_list;								/*!< The exclusion list. Atoms with the group (khist) within the list are skipped during calculations */
	};

	/*! This class is responsible for Pressure calculation */
	class Handler_Pressure : public Handler_Base {					//this class implements Handler_Base interface
	public:
		Handler_Pressure();
		virtual ~Handler_Pressure() = default;
		double get_P() const;										//!< Return the Pressure [Pa] in the system
		DPVector3 get_Pcomponents() const;							//!< Return the Pressure components {Pxx, Pyy, Pzz}
	protected:
		virtual void update_impl(System &system) override;			//!< Calculate the Pressure based on the Virial theorem

		double P;													/*!< The Pressure in the system */
		DPVector3 Pcomponents;										/*!< Pressure components {Pxx, Pyy, Pzz} */
	};

	/*! This class implements Berendsen thermostat algoritm. At each time step velocities are scaled by a scaling factor calculated based on the current and target temperatures.
		Check "Berendsen, J.Chem.Phys., 81 (1984), p. 3684" for more details. Although this method does not reproduce canonical ensemble, the average value of temperature is correct, 
		and Berendsen thermostat is widely used. To reproduce correct temperature fluctuations, Nosé–Hoover or Nosé–Hoover Chain thermostat must be used. */
	class Handler_BerendsenT : public Handler_Temperature {			//this class expands Handler_Temperature
	public:
		Handler_BerendsenT(double T, double tau = 2.0);				//!< Initialize Handler_BerendsenT object. T is the target temperature, and tau is the time constant that defines the strength of coupling to the thermal bath
		virtual ~Handler_BerendsenT() = default;
		void set_T(double T);										//!< Set the target temperature
		virtual std::string info() const override;					//!< Return a string with information about a handler
	protected:
		virtual void update_impl(System &system) override;			//!< Update particle velocities to maintain T0

		double T0;													/*!< The target temperature */
		double tau;													/*!< The time constant [ps] that defines the strength of coupling to the thermal bath. It is usually chosen from the range of 0.1 - 2 ps */
	};

	/*! This class implements Berendsen barostat algoritm. At each time step the system size is scaled to maintain the target temperature. Check "Berendsen, J.Chem.Phys., 81 (1984), p. 3684" for more details. */
	class Handler_BerendsenP : public Handler_Pressure {			//this class expands Handler_Pressure
	public:
		 /*! This is an enum class describing methods of pressure control */
		enum class Mode {XYZ,										/*!< X, Y, and Z sizes are scaled together */
			XY_Zindependent,										/*!< X and Y sizes are controlled together based on 0.5*(Pxx + Pyy), and Z size is controlled based on Pzz */
			XYZindependent,											/*!< X, Y, and Z sizes are scaled independently based on Pxx, Pyy, and Pzz */
			Z														/*!< Control only Z size based on Pzz */
		};

		Handler_BerendsenP(double P, Mode mode = Mode::XYZ, double beta = 5e-12); /**<  Initialize Handler_BerendsenP object. P [Pa] is the target temperature, mode is the method of pressure control,  
																	* and beta = Kt*dt/(3*tp) [Pa^-1] (Kt - isothermal compressibility, dt - time step, tp - characteristic time of pressure relaxation)*/
		virtual ~Handler_BerendsenP() = default;
		void set_P(double P);										//!< Set the target pressure
		virtual std::string info() const override;					//!< Return a string with information about a handler
	protected:
		virtual void update_impl(System &system) override;			//!< Update system dimensions to maintain P0
		std::string to_string(Mode mode) const;						//!< Convert Mode (pressure control method) into a string

		double P0;													/*!< The target pressure */
		double beta;												/*!< beta = Kt*dt/(3*tp) [Pa^-1] (Kt - isothermal compressibility, dt - time step, tp - characteristic time of pressure relaxation */
		Mode mode;													/*!< The method of pressure control, check enum class Mode for more details */
	};

	/*! This class is responsible for heating of the material in a less violent manner then Handler_Velocity. However, a small value of initial velocity must still be assigned */
	class Handler_Heating : public Handler_Temperature {			//this class expands Handler_Temperature
	public:
		Handler_Heating(double T, double tau = 1.0);				//!< Initialize Handler_Heating object. T is the target temperature, and tau is a characteristic time [ps] of heating
		virtual ~Handler_Heating() = default;
		virtual std::string info() const override;					//!< Return a string with information about a handler
	protected:
		virtual void update_impl(System &system) override;			//!< Perform heating of the system

		double T0;													/*!< The target temperature */
		double tau;													/*!< The characteristic time of heating [ps] */
		bool done;													/*!< A flag showing if heating is completed */
	};

	/*! This class implements a Quench algorithm that allows to stop motion of atoms and bring system to the minimum (local or global) of the potential surface */
	class Handler_Quench : public Handler_Base {					//this class implements Handler_Base interface
	public:
		Handler_Quench();
		virtual ~Handler_Quench() = default;
		virtual std::string info() const override;					//!< Return a string with information about a handler
	protected:
		virtual void update_impl(System &system) override;			//!< Quench atomic motion

		std::vector<double> Ek_prev;								/*!< The kinetic energy of each atom on the previous time step */
	};
}
