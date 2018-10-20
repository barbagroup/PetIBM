#include <petscsys.h>
#include <petscksp.h>

#include <petibm/mesh.h>
#include <petibm/type.h>

namespace petibm
{
namespace misc
{
class LinInterpBase
{
public:
    LinInterpBase() = default;

    LinInterpBase(const MPI_Comm &comm,
                  const type::RealVec1D &point,
                  const type::Mesh &mesh,
                  const type::Field &field);

    ~LinInterpBase();

    PetscErrorCode destroy();

    PetscErrorCode getValue(const DM &da, const Vec &vec, PetscReal &val);

protected:
    type::RealVec1D target;

    type::RealVec1D bl, tr;

    type::IntVec1D idxDirs;

    type::RealVec1D base_a;

    type::RealVec1D Op_a;

    Vec sub, base, coeffs;

    Mat Op;

    KSP ksp;

    virtual PetscErrorCode init(const type::RealVec1D &point,
                                const type::Mesh &mesh,
                                const type::Field &field) = 0;

    PetscErrorCode getBLGridlineIndices(const type::Mesh &mesh,
                                        const type::Field &field);

    PetscErrorCode getBoxCoords(const type::Mesh &mesh,
                                const type::Field &field);

    virtual PetscErrorCode setUpBase() = 0;

    virtual PetscErrorCode setUpOp() = 0;

    PetscErrorCode setUpKSP();

    virtual PetscErrorCode setSubVec(const DM &da, const Vec &vec) = 0;
};  // LinInterpBase

class TriLinInterp: public LinInterpBase
{
public:
    TriLinInterp(const MPI_Comm &comm,
                 const type::RealVec1D &point,
                 const type::Mesh &mesh,
                 const type::Field &field);
    ~TriLinInterp() = default;

protected:
    PetscErrorCode init(const type::RealVec1D &point,
                        const type::Mesh &mesh,
                        const type::Field &field);

    PetscErrorCode setUpBase();

    PetscErrorCode setUpOp();

    PetscErrorCode setSubVec(const DM &da, const Vec &vec);
};  // TriLinInterp

class BiLinInterp: public LinInterpBase
{
public:
    BiLinInterp(const MPI_Comm &comm,
                const type::RealVec1D &point,
                const type::Mesh &mesh,
                const type::Field &field);
    ~BiLinInterp() = default;

protected:
    PetscErrorCode init(const type::RealVec1D &point,
                        const type::Mesh &mesh,
                        const type::Field &field);

    PetscErrorCode setUpBase();

    PetscErrorCode setUpOp();

    PetscErrorCode setSubVec(const DM &da, const Vec &vec);
};  // BiLinInterp

}  // end of namespace misc

namespace type
{
typedef std::shared_ptr<misc::LinInterpBase> LinInterp;

}  // end of namespace type

namespace misc
{
PetscErrorCode createLinInterp(const MPI_Comm &comm,
                               const type::RealVec1D &point,
                               const type::Mesh &mesh,
                               const type::Field &field,
                               type::LinInterp &interp);
}  // end of namespace misc

}  // end of namespace petibm
