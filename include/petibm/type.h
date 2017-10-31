/***************************************************************************//**
 * \file type.h
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Definition of enumerated types.
 */

# pragma once

// here goes C++ STL
# include <string>
# include <vector>
# include <map>

// here goes PETSc
# include <petscsys.h>
# include <petscmat.h>


/** \brief operator < for MatStencil for using it as a key to map. 
 *
 * For the detail of implimentation, please refer to 
 * std::lexicographical_compare. This is just a splified version.
 */
bool operator<(const MatStencil &, const MatStencil &);


namespace petibm
{
/** \brief frequently used types are defined here */
namespace type
{
    /** \brief legal physical directions. */
    enum Dir { x=0, y, z };
    /** \brief mapping between `std::string` and `Dir` type */
    extern std::map<std::string, Dir> str2dir;
    /** \brief mapping between `Dir` type and `std::string` */
    extern std::map<Dir, std::string> dir2str;


    /** \brief legal fields. */
    enum Field { u=0, v, w , p, vertex};
    /** \brief mapping between `std::string` and `Field` type */
    extern std::map<std::string, Field> str2fd;
    /** \brief mapping between `Field` type and `std::string` */
    extern std::map<Field, std::string> fd2str;


    /** \brief type of boundary condition. */
    enum BCType { NOBC=0, PERIODIC, DIRICHLET, NEUMANN, CONVECTIVE };
    /** \brief mapping between `std::string` and `BCType` type */
    extern std::map<std::string, BCType> str2bt;
    /** \brief mapping between `BCType` type and `std::string` */
    extern std::map<BCType, std::string> bt2str;


    /** \brief location of the boundary. */
    enum BCLoc { XMINUS=0, XPLUS, YMINUS, YPLUS, ZMINUS, ZPLUS };
    /** \brief mapping between `std::string` and `BCLoc` type */
    extern std::map<std::string, BCLoc> str2bl;
    /** \brief mapping between `BCLoc` type and `std::string` */
    extern std::map<BCLoc, std::string> bl2str;


    /** \brief Immersed boundary method used to solve the flow. */
    enum IBMethod
    {
        NAVIER_STOKES = 0,  ///< no immersed boundary (Perot, 1993)
        TAIRA_COLONIUS = 1  ///< immersed boundary projection method (Taira & Colonius, 2007)
    };
    /** \brief mapping between `std::string` and `IBMethod` type */
    extern std::map<std::string, IBMethod> str2ibm;
    /** \brief mapping between `IBMethod` type and `std::string` */
    extern std::map<IBMethod, std::string> ibm2str;


    /** \brief Staggered mode to define the location of mesh points. */
    enum StaggeredMode
    {
        CELL_CENTERED = 0,     ///< cell-centered quantity
        STAGGERED_MODE_X = 1,  ///< x-component of staggered quantity
        STAGGERED_MODE_Y = 2,  ///< y-component of staggered quantity
        STAGGERED_MODE_Z = 3   ///< z-component of staggered quantity
    };
    /** \brief mapping between `std::string` and `StaggeredMode` type */
    extern std::map<std::string, StaggeredMode> str2sm;
    /** \brief mapping between `StaggeredMode` type and `std::string` */
    extern std::map<StaggeredMode, std::string> sm2str;


    /** \brief output format */
    enum OutputType {Binary = 0, VTK = 1, HDF5 = 2};
    /** \brief mapping between `std::string` and `OutputType` type */
    extern std::map<std::string, OutputType> str2out;
    /** \brief mapping between `OutputType` type and `std::string` */
    extern std::map<OutputType, std::string> out2str;


    /** \brief 1D std::vector holding PetscInt. */
    typedef std::vector<PetscInt>   IntVec1D;
    /** \brief 2D std::vector holding PetscInt. */
    typedef std::vector<IntVec1D>   IntVec2D;
    /** \brief 3D std::vector holding PetscInt. */
    typedef std::vector<IntVec2D>   IntVec3D;


    /** \brief 1D std::vector holding PetscReal. */
    typedef std::vector<PetscReal>  RealVec1D;
    /** \brief 2D std::vector holding PetscReal. */
    typedef std::vector<RealVec1D>  RealVec2D;
    /** \brief 3D std::vector holding PetscReal. */
    typedef std::vector<RealVec2D>  RealVec3D;


    /** \brief 1D std::vector holding PetscBool. */
    typedef std::vector<PetscBool>  BoolVec1D;
    /** \brief 2D std::vector holding PetscBool. */
    typedef std::vector<BoolVec1D>  BoolVec2D;
    /** \brief 3D std::vector holding PetscBool. */
    typedef std::vector<BoolVec2D>  BoolVec3D;


    /** \brief a vector of pointer to mimic a collection of ghosted 1D vectors. */
    typedef std::vector<PetscReal*>             GhostedVec2D;
    /** \brief a vector of vector pointer to mimic a collection of ghosted 2D vectors. */
    typedef std::vector<GhostedVec2D>           GhostedVec3D;


    /** \brief a structure holding information of perturbation. */
    struct Perturbation
    {
        PetscReal   freq = 0.;  ///< frequency
        PetscReal   amp = 0.;  ///< amplitude
    };


    /** \brief a struct holding information of output */
    struct OutputInfo
    {
        OutputType      format = OutputType::Binary;
        PetscBool       outputFlux = PETSC_TRUE;
        PetscBool       outputVelocity = PETSC_FALSE;
    };


    /** \brief a struct holding information of time-stepping */
    struct SteppingInfo
    {
        PetscReal       dt;
        PetscInt        nStart;
        PetscInt        nTotal;
        PetscInt        nSave;
        PetscInt        nRestart;
    };


    /** \brief a data structure for a single ghost point. */
    struct GhostPointInfo
    {
        PetscInt    lclId;
        MatStencil  targetStencil;
        PetscInt    targetPackedId;
        PetscReal   area, dL;
        PetscReal   a0, a1;
        PetscReal   value;
    };


    /** \brief a map (dictionary) for ghost points. */
    typedef std::map<MatStencil, GhostPointInfo> GhostPointsList;


    /** \brief a struct holding information about which row in a matrix should
     *         be modified based on BCs.
     */
    struct RowModifier
    {
        PetscInt    row;
        PetscReal   coeff;
    };


    /** \brief a type that holds necessary for a matrix modifier that modifies
     *         matrix coefficient based on BCs.
     */
    typedef std::vector<std::map<MatStencil, RowModifier>> MatrixModifier;

} // end of namespace type
} // end of namespace petibm

/**
 * \class PetscInt
 * \brief PETSc type that represents an integer.
 * \sa http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscInt.html
 */

/**
 * \class PetscReal
 * \brief PETSc type that represents a real number.
 * \sa http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscReal.html
 */

/**
 * \class PetscBool
 * \brief PETSc type that represents a logical variable.
 * \sa http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscBool.html
 */

/**
 * \class PetscErrorCode
 * \brief PETSc type used to return the error code from a function.
 * \sa http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscErrorCode.html
 */

/**
 * \var EULER_EXPLICIT
 * \f[ \frac{u^{n+1}-u^n}{\Delta t} = f(u^n) \f]
 */

/**
 * \var EULER_IMPLICIT
 * \f[ \frac{u^{n+1}-u^n}{\Delta t} = f(u^{n+1}) \f]
 */

/**
 * \var ADAMS_BASHFORTH_2
 * \f[ \frac{u^{n+1}-u^n}{\Delta t} = \frac{3}{2}f(u^n) - \frac{1}{2}f(u^{n-1}) \f]
 */

/**
 * \var CRANK_NICOLSON
 * \f[ \frac{u^{n+1}-u^n}{\Delta t} = \frac{1}{2}\left(f(u^{n+1}) + f(u^n)\right) \f]
 */

/**
 * \var NAVIER_STOKES
 * Solve the incompressible Navier-Stokes equations in the domain without the
 * presence of any immersed boundaries. The method used is the fractional step
 * method as formulated in the paper by Blair Perot
 * <a href="http://dx.doi.org/10.1006/jcph.1993.1162">(J. Comput. Phys, 108.1, 1993)</a>
 */

/**
 * \var TAIRA_COLONIUS
 * Solve the flow using the immersed boundary projection method proposed by
 * Taira and Colonius
 * <a href="http://dx.doi.org/10.1016/j.jcp.2007.03.005">(J. Comput. Phys, 225.2, 2007)</a>
 */