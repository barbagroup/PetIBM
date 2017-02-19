/***************************************************************************//**
 * \file types.cpp
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Implementations of members inside namespace `types`.
 */

// here goes our headers
#include "types.h"


/** \copydoc types */
namespace types
{
    /** \copydoc types::str2dir */
    std::map<std::string, Dir> str2dir {{"x", x}, {"y", y}, {"z", z}};
    /** \copydoc types::dir2str */
    std::map<Dir, std::string> dir2str {{x, "x"}, {y, "y"}, {z, "z"}};


    /** \copydoc types::str2vc */
    std::map<std::string, VelocityComponent> str2vc {{"u", u}, {"v", v}, {"w", w}};
    /** \copydoc types::vc2str */
    std::map<VelocityComponent, std::string> vc2str {{u, "u"}, {v, "v"}, {w, "w"}};


    /** \copydoc types::str2bt */
    std::map<std::string, BCType> str2bt {
        {"DIRICHLET", DIRICHLET}, {"NEUMANN", NEUMANN}, 
        {"CONVECTIVE", CONVECTIVE}, {"PERIODIC", PERIODIC}};
    /** \copydoc types::bt2str */
    std::map<BCType, std::string> bt2str {
        {DIRICHLET, "DIRICHLET"}, {NEUMANN, "NEUMANN"}, 
        {CONVECTIVE, "CONVECTIVE"}, {PERIODIC, "PERIODIC"}};

    /** \copydoc types::str2bl */
    std::map<std::string, BCLoc> str2bl {
        {"left", XMINUS}, {"right", XPLUS}, {"bottom", YMINUS}, 
        {"top", YPLUS}, {"back", ZMINUS}, {"front", ZPLUS},
        {"xMinus", XMINUS}, {"xPlus", XPLUS}, {"yMinus", YMINUS}, 
        {"yPlus", YPLUS}, {"zMinus", ZMINUS}, {"zPlus", ZPLUS},
        {"XMINUS", XMINUS}, {"XPLUS", XPLUS}, {"YMINUS", YMINUS}, 
        {"YPLUS", YPLUS}, {"ZMINUS", ZMINUS}, {"ZPLUS", ZPLUS}};
    /** \copydoc types::bl2str */
    std::map<BCLoc, std::string> bl2str {
        {XMINUS, "xMinus"}, {XPLUS, "xPlus"}, {YMINUS, "yMinus"}, 
        {YPLUS, "yPlus"}, {ZMINUS, "zMinus"}, {ZPLUS, "zPlus"}};


    /** \copydoc types::str2ts */
    std::map<std::string, TimeScheme> str2ts { {"NONE", NONE}, 
        {"EULER_EXPLICIT", EULER_EXPLICIT}, {"EULER_IMPLICIT", EULER_IMPLICIT}, 
        {"ADAMS_BASHFORTH_2", ADAMS_BASHFORTH_2}, {"CRANK_NICOLSON", CRANK_NICOLSON}};
    /** \copydoc types::ts2str */
    std::map<TimeScheme, std::string> ts2str { {NONE, "NONE"}, 
        {EULER_EXPLICIT, "EULER_EXPLICIT"}, {EULER_IMPLICIT, "EULER_IMPLICIT"}, 
        {ADAMS_BASHFORTH_2, "ADAMS_BASHFORTH_2"}, {CRANK_NICOLSON, "CRANK_NICOLSON"}};


    /** \copydoc types::str2ibm */
    std::map<std::string, IBMethod> str2ibm {
        {"NAVIER_STOKES", NAVIER_STOKES}, {"TAIRA_COLONIUS", TAIRA_COLONIUS}};
    /** \copydoc types::ibm2str */
    std::map<IBMethod, std::string> ibm2str {
        {NAVIER_STOKES, "NAVIER_STOKES"}, {TAIRA_COLONIUS, "TAIRA_COLONIUS"}};


    /** \copydoc types::str2sm */
    std::map<std::string, StaggeredMode> str2sm {
        {"CELL_CENTERED", CELL_CENTERED}, {"STAGGERED_MODE_X", STAGGERED_MODE_X}, 
        {"STAGGERED_MODE_Y", STAGGERED_MODE_Y}, {"STAGGERED_MODE_Z", STAGGERED_MODE_Z}};
    /** \copydoc types::sm2str */
    std::map<StaggeredMode, std::string> sm2str {
        {CELL_CENTERED, "CELL_CENTERED"}, {STAGGERED_MODE_X, "STAGGERED_MODE_X"}, 
        {STAGGERED_MODE_Y, "STAGGERED_MODE_Y"}, {STAGGERED_MODE_Z, "STAGGERED_MODE_Z"}};


    /** \copydoc types::str2et */
    std::map<std::string, ExecuteType> str2et {{"CPU", CPU}, {"GPU", GPU}};
    /** \copydoc types::et2str */
    std::map<ExecuteType, std::string> et2str {{CPU, "CPU"}, {GPU, "GPU"}};
}
