/***************************************************************************//**
 * \file type.h
 * \brief Definition of user-defined types for convenience.
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \copyright MIT.
 */

# pragma once

// here goes C++ STL
# include <string>
# include <vector>
# include <map>

// here goes PETSc
# include <petscsys.h>
# include <petscmat.h>


/** 
 * \defgroup petibm PetIBM building blocks 
 * 
 * This part contains building blocks that are necessary for a flow solver in 
 * Perot (1993) framework and staggered-grid finite difference method.
 */


/** 
 * \defgroup type User-defined types in PetIBM
 * \brief Useful user-defined types in PetIBM.
 * 
 * User-defined types to simplify the source code.
 * 
 * PETSc also has its own type definitions. Here are information of some PETSc
 * types we use in PetIBM:
 * 
 * \arg \c \b PetscInt http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscInt.html
 * \arg \c \b PetscReal http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscReal.html
 * \arg \c \b PetscBool http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscBool.html
 * \arg \c \b PetscErrorCode http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscErrorCode.html
 * \arg \c \b MatStencil http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatStencil.html
 * \arg \c \b Mat http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/Mat.html
 * \arg \c \b Vec http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/Vec.html
 * \arg \c \b KSP http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSP.html
 * \ingroup petibm
 */


/** 
 * \brief Operator< of MatStencil for using it as a key to map. 
 *
 * For the detail of implementation, please refer to 
 * std::lexicographical_compare. This is just a simplified version.
 */
bool operator<(const MatStencil &, const MatStencil &);


/** 
 * \brief A toolbox for building flow solvers. 
 * \ingroup petibm
 */
namespace petibm
{
    
/** 
 * \brief Frequently used types, structures, and enums.
 * \ingroup type
 */
namespace type
{
    /** \brief Legal physical directions. \ingroup type */
    enum Dir { x=0, y, z };
    /** \brief Mapping between `std::string` and \ref Dir. \ingroup type */
    extern std::map<std::string, Dir> str2dir;
    /** \brief Mapping between \ref Dir and `std::string`. \ingroup type */
    extern std::map<Dir, std::string> dir2str;


    /** \brief Legal fields. \ingroup type */
    enum Field { u=0, v, w, p, vertex};
    /** \brief Mapping between `std::string` and \ref Field. \ingroup type */
    extern std::map<std::string, Field> str2fd;
    /** \brief Mapping between \ref Field and `std::string`. \ingroup type */
    extern std::map<Field, std::string> fd2str;


    /** \brief Type of boundary conditions. \ingroup type */
    enum BCType { NOBC=0, PERIODIC, DIRICHLET, NEUMANN, CONVECTIVE };
    /** \brief Mapping between `std::string` and \ref BCType. \ingroup type */
    extern std::map<std::string, BCType> str2bt;
    /** \brief mapping between \ref BCType and `std::string`. \ingroup type */
    extern std::map<BCType, std::string> bt2str;


    /** \brief Location of a boundary. \ingroup type */
    enum BCLoc { XMINUS=0, XPLUS, YMINUS, YPLUS, ZMINUS, ZPLUS };
    /** \brief Mapping between `std::string` and \ref BCLoc. \ingroup type */
    extern std::map<std::string, BCLoc> str2bl;
    /** \brief Mapping between \ref BCLoc and `std::string`. \ingroup type */
    extern std::map<BCLoc, std::string> bl2str;


    /** \brief 1D std::vector holding PetscInt. \ingroup type */
    typedef std::vector<PetscInt>   IntVec1D;
    /** \brief 2D std::vector holding PetscInt. \ingroup type */
    typedef std::vector<IntVec1D>   IntVec2D;
    /** \brief 3D std::vector holding PetscInt. \ingroup type */
    typedef std::vector<IntVec2D>   IntVec3D;


    /** \brief 1D std::vector holding PetscReal. \ingroup type */
    typedef std::vector<PetscReal>  RealVec1D;
    /** \brief 2D std::vector holding PetscReal. \ingroup type */
    typedef std::vector<RealVec1D>  RealVec2D;
    /** \brief 3D std::vector holding PetscReal. \ingroup type */
    typedef std::vector<RealVec2D>  RealVec3D;


    /** \brief 1D std::vector holding PetscBool. \ingroup type */
    typedef std::vector<PetscBool>  BoolVec1D;
    /** \brief 2D std::vector holding PetscBool. \ingroup type */
    typedef std::vector<BoolVec1D>  BoolVec2D;
    /** \brief 3D std::vector holding PetscBool. \ingroup type */
    typedef std::vector<BoolVec2D>  BoolVec3D;


    /** \brief a vector of pointers to mimic ghosted 1D vectors. \ingroup type */
    typedef std::vector<PetscReal*>             GhostedVec2D;
    /** \brief a vector of vector pointers to mimic ghosted 2D vectors. \ingroup type */
    typedef std::vector<GhostedVec2D>           GhostedVec3D;


    /** \brief A data structure for a single ghost point. \ingroup type */
    struct GhostPointInfo
    {
        PetscInt    lclId; ///< the index in a local velocity vector of this ghost point
        MatStencil  targetStencil; ///< the stencil of point corresponding to this ghost point
        PetscInt    targetPackedId; ///< the index of target point in global packed velocity vector
        PetscReal   area; ///< the flux area of this ghost point
        PetscReal   dL; ///< the distance between this ghost point and target velocity point
        PetscReal   a0; ///< the coefficient of the velocity at the target point
        PetscReal   a1; ///< the constant coefficient
        PetscReal   value; ///< the value of this ghost point
    };


    /** \brief A map between MatStencil and GhostPointInfo. \ingroup type */
    typedef std::map<MatStencil, GhostPointInfo> GhostPointsList;


    /** 
     * \brief A struct holding information about which row in a matrix should 
     * be modified based on BCs.
     * \ingroup type 
     */
    struct RowModifier
    {
        PetscInt    row; ///< the global index of a row in a Mat
        PetscReal   coeff; ///< coefficient
    };


    /**
     * \brief A type that holds necessary info for a matrix modifier that modifies
     * matrix coefficient based on BCs.
     * \ingroup type
     */
    typedef std::vector<std::map<MatStencil, RowModifier>> MatrixModifier;

} // end of namespace type
} // end of namespace petibm
