/**
 * @file StFlowLe.h
 */

#ifndef CT_STFLOWLE_H
#define CT_STFLOWLE_H
#include "cantera/transport.h"       // XXX ars
#include "cantera/numerics/funcs.h"  // XXX ars
#include "cantera/base/ctml.h"       // XXX ars
//#include "cantera/oneD/Domain1D.h"
#include "Domain1DLe.h"
#include "cantera/base/Array.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/kinetics/Kinetics.h"

// XXX ars
//using namespace Cantera;
using Cantera::IdealGasPhase;
using Cantera::thermo_t;
using Cantera::Kinetics;
using Cantera::Transport;
using Cantera::Array2D;
using Cantera::warn_deprecated;
using Cantera::Pi;
using Cantera::GasConstant;
#ifdef OLD_CANTERA
using Cantera::cMulticomponent;
using Cantera::CK_Multicomponent;
using Cantera::cMixtureAveraged;
using Cantera::CK_MixtureAveraged;
#endif
using Cantera::linearInterp;
using Cantera::OneAtm;
using Cantera::StefanBoltz;
using Cantera::writelog;
using Cantera::writeline;
using Cantera::getFloat;
using Cantera::getFloatArray;
using Cantera::Undef;

namespace CanteraLe
{
//------------------------------------------
//   constants
//------------------------------------------

 // XXX ars: add _Le
// Offsets of solution components in the solution array.
const size_t c_offset_U_Le = 0;    // axial velocity
const size_t c_offset_V_Le = 1;    // strain rate
const size_t c_offset_T_Le = 2;    // temperature
const size_t c_offset_L_Le = 3;    // (1/r)dP/dr
const size_t c_offset_Y_Le = 4;    // mass fractions
// Transport option flags
const int c_Mixav_Transport_Le = 0;
const int c_Multi_Transport_Le = 1;
const int c_Soret_Le = 2;

//vector_fp m_Le; // XXX ars
//double* m_Le; // XXX ars

//class Transport;

/**
 *  This class represents 1D flow domains that satisfy the one-dimensional
 *  similarity solution for chemically-reacting, axisymmetric, flows.
 *  @ingroup onedim
 */
class StFlowLe : public Domain1D
{
public:
    //--------------------------------
    // construction and destruction
    //--------------------------------

    //! Create a new flow domain.
    //! @param ph Object representing the gas phase. This object will be used
    //!     to evaluate all thermodynamic, kinetic, and transport properties.
    //! @param nsp Number of species.
    //! @param points Initial number of grid points
    StFlowLe(IdealGasPhase* ph = 0, size_t nsp = 1, size_t points = 1);

    //! @name Problem Specification
    //! @{

    virtual void setupGrid(size_t n, const doublereal* z);

    thermo_t& phase() {
        return *m_thermo;
    }
    Kinetics& kinetics() {
        return *m_kin;
    }

    virtual void init() {
    }

    /**
     * Set the thermo manager. Note that the flow equations assume
     * the ideal gas equation.
     */
    void setThermo(IdealGasPhase& th) {
        m_thermo = &th;
    }

    //! Set the kinetics manager. The kinetics manager must
    void setKinetics(Kinetics& kin) {
        m_kin = &kin;
    }

    // XXX ars
    //! set the transport manager
    void setTransport(Transport& trans, const double *Le, bool withSoret = false);
    void enableSoret(bool withSoret);
    bool withSoret() const {
        return m_do_soret;
    }

    //! Set the pressure. Since the flow equations are for the limit of
    //! small Mach number, the pressure is very nearly constant
    //! throughout the flow.
    void setPressure(doublereal p) {
        m_press = p;
    }

    //! The current pressure [Pa].
    doublereal pressure() const {
        return m_press;
    }

    //! Write the initial solution estimate into array x.
    virtual void _getInitialSoln(doublereal* x) {
        for (size_t j = 0; j < m_points; j++) {
            T(x,j) = m_thermo->temperature();
            m_thermo->getMassFractions(&Y(x, 0, j));
        }
    }

    virtual void _finalize(const doublereal* x);

    //! Sometimes it is desired to carry out the simulation using a specified
    //! temperature profile, rather than computing it by solving the energy
    //! equation. This method specifies this profile.
    void setFixedTempProfile(vector_fp& zfixed, vector_fp& tfixed) {
        m_zfix = zfixed;
        m_tfix = tfixed;
    }

    /*!
     * Set the temperature fixed point at grid point j, and disable the energy
     * equation so that the solution will be held to this value.
     */
    void setTemperature(size_t j, doublereal t) {
        m_fixedtemp[j] = t;
        m_do_energy[j] = false;
    }

    /*!
     * Set the mass fraction fixed point for species k at grid point j, and
     * disable the species equation so that the solution will be held to this
     * value. Note: in practice, the species are hardly ever held fixed.
     */
    void setMassFraction(size_t j, size_t k, doublereal y) {
        m_fixedy(k,j) = y;
        m_do_species[k] = true;
    }

    //! The fixed temperature value at point j.
    doublereal T_fixed(size_t j) const {
        return m_fixedtemp[j];
    }

    //! The fixed mass fraction value of species k at point j.
    //! @deprecated Unused. To be removed after Cantera 2.2.
    doublereal Y_fixed(size_t k, size_t j) const {
        warn_deprecated("StFlowLe::Y_fixed", "To be removed after Cantera 2.2.");
        return m_fixedy(k,j);
    }

    // @}

    virtual std::string componentName(size_t n) const;

    virtual size_t componentIndex(const std::string& name) const;

    //! Print the solution.
    virtual void showSolution(const doublereal* x);

    //! Save the current solution for this domain into an XML_Node
    /*!
     *  @param o    XML_Node to save the solution to.
     *  @param sol  Current value of the solution vector. The object will pick
     *              out which part of the solution vector pertains to this
     *              object.
     */
    virtual XML_Node& save(XML_Node& o, const doublereal* const sol);

    virtual void restore(const XML_Node& dom, doublereal* soln,
                         int loglevel);

    // overloaded in subclasses
    virtual std::string flowType() {
        return "<none>";
    }

    void solveEnergyEqn(size_t j=npos) {
        bool changed = false;
        if (j == npos)
            for (size_t i = 0; i < m_points; i++) {
                if (!m_do_energy[i]) {
                    changed = true;
                }
                m_do_energy[i] = true;
            }
        else {
            if (!m_do_energy[j]) {
                changed = true;
            }
            m_do_energy[j] = true;
        }
        m_refiner->setActive(0, true);
        m_refiner->setActive(1, true);
        m_refiner->setActive(2, true);
        if (changed) {
            needJacUpdate();
        }
    }

    //! Turn radiation on / off.
    /*!
     *  The simple radiation model used was established by Y. Liu and B. Rogg
     *  [Y. Liu and B. Rogg, Modelling of thermally radiating diffusion flames
     *  with detailed chemistry and transport, EUROTHERM Seminars, 17:114-127,
     *  1991]. This model considers the radiation of CO2 and H2O.
     */
    void enableRadiation(bool doRadiation) {
        m_do_radiation = doRadiation;
    }

    //! Returns `true` if the radiation term in the energy equation is enabled
    bool radiationEnabled() const {
        return m_do_radiation;
    }

    //! Set the emissivities for the boundary values
    /*!
    *   Reads the emissivities for the left and right boundary values in the
    *   radiative term and writes them into the variables, which are used for
    *   the calculation.
    */
    void setBoundaryEmissivities(doublereal e_left, doublereal e_right) {
        if (e_left < 0 || e_left > 1) {
            throw CanteraError("setBoundaryEmissivities",
                "The left boundary emissivity must be between 0.0 and 1.0!");
        } else if (e_right < 0 || e_right > 1) {
            throw CanteraError("setBoundaryEmissivities",
                "The right boundary emissivity must be between 0.0 and 1.0!");
        } else {
            m_epsilon_left = e_left;
            m_epsilon_right = e_right;
        }
    }

    void fixTemperature(size_t j=npos) {
        bool changed = false;
        if (j == npos)
            for (size_t i = 0; i < m_points; i++) {
                if (m_do_energy[i]) {
                    changed = true;
                }
                m_do_energy[i] = false;
            }
        else {
            if (m_do_energy[j]) {
                changed = true;
            }
            m_do_energy[j] = false;
        }
        m_refiner->setActive(0, false);
        m_refiner->setActive(1, false);
        m_refiner->setActive(2, false);
        if (changed) {
            needJacUpdate();
        }
    }

    //! @deprecated Species equations are always solved. To be removed after
    //! Cantera 2.2.
    bool doSpecies(size_t k) {
        warn_deprecated("StFlowLe::doSpecies", "To be removed after Cantera 2.2.");
        return m_do_species[k];
    }
    bool doEnergy(size_t j) {
        return m_do_energy[j];
    }

    //! @deprecated Species equations are always solved. To be removed after
    //! Cantera 2.2.
    void solveSpecies(size_t k=npos) {
        warn_deprecated("StFlowLe::solveSpecies", "To be removed after Cantera 2.2.");
        if (k == npos) {
            for (size_t i = 0; i < m_nsp; i++) {
                m_do_species[i] = true;
            }
        } else {
            m_do_species[k] = true;
        }
        needJacUpdate();
    }

    //! @deprecated Species equations are always solved. To be removed after
    //! Cantera 2.2.
    void fixSpecies(size_t k=npos) {
        warn_deprecated("StFlowLe::fixSpecies", "To be removed after Cantera 2.2.");
        if (k == npos) {
            for (size_t i = 0; i < m_nsp; i++) {
                m_do_species[i] = false;
            }
        } else {
            m_do_species[k] = false;
        }
        needJacUpdate();
    }

    void integrateChem(doublereal* x,doublereal dt);

    //! Change the grid size. Called after grid refinement.
    void resize(size_t components, size_t points);

    virtual void setFixedPoint(int j0, doublereal t0) {}

    void setJac(MultiJac* jac);

    //! Set the gas object state to be consistent with the solution at point j.
    void setGas(const doublereal* x, size_t j);

    //! Set the gas state to be consistent with the solution at the midpoint
    //! between j and j + 1.
    void setGasAtMidpoint(const doublereal* x, size_t j);

    doublereal density(size_t j) const {
        return m_rho[j];
    }

    virtual bool fixed_mdot() {
        return true;
    }
    void setViscosityFlag(bool dovisc) {
        m_dovisc = dovisc;
    }

    /*!
     *  Evaluate the residual function for axisymmetric stagnation flow. If
     *  jpt is less than zero, the residual function is evaluated at all grid
     *  points. If jpt >= 0, then the residual function is only evaluated at
     *  grid points jpt-1, jpt, and jpt+1. This option is used to efficiently
     *  evaluate the Jacobian numerically.
     */
    virtual void eval(size_t j, doublereal* x, doublereal* r,
                      integer* mask, doublereal rdt);

    //! Evaluate all residual components at the right boundary.
    virtual void evalRightBoundary(doublereal* x, doublereal* res,
                                   integer* diag, doublereal rdt) = 0;

    //! Evaluate the residual corresponding to the continuity equation at all
    //! interior grid points.
    virtual void evalContinuity(size_t j, doublereal* x, doublereal* r,
                                integer* diag, doublereal rdt) = 0;

protected:
    doublereal component(const doublereal* x, size_t i, size_t j) const {
        return x[index(i,j)];
    }

    doublereal conc(const doublereal* x, size_t k,size_t j) const {
        return Y(x,k,j)*density(j)/m_wt[k];
    }

    doublereal cbar(const doublereal* x, size_t k, size_t j) const {
        return std::sqrt(8.0*GasConstant * T(x,j) / (Pi * m_wt[k]));
    }

    doublereal wdot(size_t k, size_t j) const {
        return m_wdot(k,j);
    }

    //! Write the net production rates at point `j` into array `m_wdot`
    void getWdot(doublereal* x, size_t j) {
        setGas(x,j);
        m_kin->getNetProductionRates(&m_wdot(0,j));
    }

    /**
     * Update the thermodynamic properties from point j0 to point j1
     * (inclusive), based on solution x.
     */
    void updateThermo(const doublereal* x, size_t j0, size_t j1) {
        for (size_t j = j0; j <= j1; j++) {
            setGas(x,j);
            m_rho[j] = m_thermo->density();
            m_wtm[j] = m_thermo->meanMolecularWeight();
            m_cp[j]  = m_thermo->cp_mass();
        }
    }

    //--------------------------------
    // central-differenced derivatives
    //--------------------------------

    doublereal cdif2(const doublereal* x, size_t n, size_t j,
                     const doublereal* f) const {
        doublereal c1 = (f[j] + f[j-1])*(x[index(n,j)] - x[index(n,j-1)]);
        doublereal c2 = (f[j+1] + f[j])*(x[index(n,j+1)] - x[index(n,j)]);
        return (c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
    }


    //! @name Solution components
    //! @{

    doublereal T(const doublereal* x, size_t j) const {
        return x[index(c_offset_T_Le, j)];
    }
    doublereal& T(doublereal* x, size_t j) {
        return x[index(c_offset_T_Le, j)];
    }
    doublereal T_prev(size_t j) const {
        return prevSoln(c_offset_T_Le, j);
    }

    doublereal rho_u(const doublereal* x, size_t j) const {
        return m_rho[j]*x[index(c_offset_U_Le, j)];
    }

    doublereal u(const doublereal* x, size_t j) const {
        return x[index(c_offset_U_Le, j)];
    }

    doublereal V(const doublereal* x, size_t j) const {
        return x[index(c_offset_V_Le, j)];
    }
    doublereal V_prev(size_t j) const {
        return prevSoln(c_offset_V_Le, j);
    }

    doublereal lambda(const doublereal* x, size_t j) const {
        return x[index(c_offset_L_Le, j)];
    }

    doublereal Y(const doublereal* x, size_t k, size_t j) const {
        return x[index(c_offset_Y_Le + k, j)];
    }

    doublereal& Y(doublereal* x, size_t k, size_t j) {
        return x[index(c_offset_Y_Le + k, j)];
    }

    doublereal Y_prev(size_t k, size_t j) const {
        return prevSoln(c_offset_Y_Le + k, j);
    }

    doublereal X(const doublereal* x, size_t k, size_t j) const {
        return m_wtm[j]*Y(x,k,j)/m_wt[k];
    }

    doublereal flux(size_t k, size_t j) const {
        return m_flux(k, j);
    }
    //! @}

    //! @name convective spatial derivatives.
    //! These use upwind differencing, assuming u(z) is negative
    //! @{
    doublereal dVdz(const doublereal* x, size_t j) const {
        size_t jloc = (u(x,j) > 0.0 ? j : j + 1);
        return (V(x,jloc) - V(x,jloc-1))/m_dz[jloc-1];
    }

    doublereal dYdz(const doublereal* x, size_t k, size_t j) const {
        size_t jloc = (u(x,j) > 0.0 ? j : j + 1);
        return (Y(x,k,jloc) - Y(x,k,jloc-1))/m_dz[jloc-1];
    }

    doublereal dTdz(const doublereal* x, size_t j) const {
        size_t jloc = (u(x,j) > 0.0 ? j : j + 1);
        return (T(x,jloc) - T(x,jloc-1))/m_dz[jloc-1];
    }
    //! @}

    doublereal shear(const doublereal* x, size_t j) const {
        doublereal c1 = m_visc[j-1]*(V(x,j) - V(x,j-1));
        doublereal c2 = m_visc[j]*(V(x,j+1) - V(x,j));
        return 2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
    }

    doublereal divHeatFlux(const doublereal* x, size_t j) const {
        doublereal c1 = m_tcon[j-1]*(T(x,j) - T(x,j-1));
        doublereal c2 = m_tcon[j]*(T(x,j+1) - T(x,j));
        return -2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
    }

    size_t mindex(size_t k, size_t j, size_t m) {
        return m*m_nsp*m_nsp + m_nsp*j + k;
    }

    //! Update the diffusive mass fluxes.
    void updateDiffFluxes(const doublereal* x, size_t j0, size_t j1);

    //---------------------------------------------------------
    //             member data
    //---------------------------------------------------------

    doublereal m_press;        // pressure

    // grid parameters
    vector_fp m_dz;

    // mixture thermo properties
    vector_fp m_rho;
    vector_fp m_wtm;

    // species thermo properties
    vector_fp m_wt;
    vector_fp m_cp;

    // transport properties
    vector_fp m_visc;
    vector_fp m_tcon;
    vector_fp m_diff;
    vector_fp m_multidiff;
    Array2D m_dthermal;
    Array2D m_flux;
    vector_fp m_Le; // XXX ars

    // production rates
    Array2D m_wdot;

    size_t m_nsp;

    IdealGasPhase* m_thermo;
    Kinetics* m_kin;
    Transport* m_trans;

    MultiJac* m_jac;

    // boundary emissivities for the radiation calculations
    doublereal m_epsilon_left;
    doublereal m_epsilon_right;

    //! Indices within the ThermoPhase of the radiating species. First index is
    //! for CO2, second is for H2O.
    std::vector<size_t> m_kRadiating;

    // flags
    std::vector<bool> m_do_energy;
    bool m_do_soret;
    std::vector<bool> m_do_species;
    int m_transport_option;

    // flag for the radiative heat loss
    bool m_do_radiation;

    // radiative heat loss vector
    // vector which contains the values of the radiative heat loss
    vector_fp m_qdotRadiation;

    // fixed T and Y values
    Array2D   m_fixedy; //!< @deprecated
    vector_fp m_fixedtemp;
    vector_fp m_zfix;
    vector_fp m_tfix;

    bool m_dovisc;

    //! Update the transport properties at grid points in the range from `j0`
    //! to `j1`, based on solution `x`.
    void updateTransport(doublereal* x, size_t j0, size_t j1);

private:
    vector_fp m_ybar;
};

/**
 * A class for axisymmetric stagnation flows.
 * @ingroup onedim
 */
class AxiStagnFlowLe : public StFlowLe
{
public:
    AxiStagnFlowLe(IdealGasPhase* ph = 0, size_t nsp = 1, size_t points = 1) :
        StFlowLe(ph, nsp, points) {
        m_dovisc = true;
    }

    virtual void evalRightBoundary(doublereal* x, doublereal* res,
                                   integer* diag, doublereal rdt);
    virtual void evalContinuity(size_t j, doublereal* x, doublereal* r,
                                integer* diag, doublereal rdt);

    virtual std::string flowType() {
        return "Axisymmetric Stagnation";
    }
};

/**
 * A class for freely-propagating premixed flames.
 * @ingroup onedim
 */
class FreeFlameLe : public StFlowLe
{
public:
    FreeFlameLe(IdealGasPhase* ph = 0, size_t nsp = 1, size_t points = 1);
    virtual void evalRightBoundary(doublereal* x, doublereal* res,
                                   integer* diag, doublereal rdt);
    virtual void evalContinuity(size_t j, doublereal* x, doublereal* r,
                                integer* diag, doublereal rdt);

    virtual std::string flowType() {
        return "Free Flame";
    }
    virtual bool fixed_mdot() {
        return false;
    }
    virtual void _finalize(const doublereal* x);
    virtual void restore(const XML_Node& dom, doublereal* soln, int loglevel);

    virtual XML_Node& save(XML_Node& o, const doublereal* const sol);

    //! Location of the point where temperature is fixed
    doublereal m_zfixed;

    //! Temperature at the point used to fix the flame location
    doublereal m_tfixed;
};

}

#endif
