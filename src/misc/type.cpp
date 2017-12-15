/***************************************************************************//**
 * \file type.cpp
 * \author Anush Krishnan (anush@bu.edu)
 * \author Olivier Mesnard (mesnardo@gwu.edu)
 * \author Pi-Yueh Chuang (pychuang@gwu.edu)
 * \brief Implementations of members inside namespace \ref petibm::type.
 * \copyright MIT.
 */

// here goes our headers
#include <petibm/type.h>


// implement operator< for MatStencil
bool operator<(const MatStencil &l, const MatStencil &r)
{
    if (l.k < r.k) return true;
    if (l.k > r.k) return false;

    if (l.j < r.j) return true;
    if (l.j > r.j) return false;

    if (l.i < r.i) return true;
    if (l.i > r.i) return false;

    if (l.c < r.c) return true;

    return false;
}


namespace petibm
{
namespace type
{
    // default values of type::str2dir
    std::map<std::string, Dir> str2dir {{"x", x}, {"y", y}, {"z", z}};
    // default values of type::dir2str
    std::map<Dir, std::string> dir2str {{x, "x"}, {y, "y"}, {z, "z"}};

    // default values of type::str2fd
    std::map<std::string, Field> str2fd {
        {"u", u}, {"v", v}, {"w", w}, {"p", p}, {"vertex", vertex}};
    // default values of type::fd2str
    std::map<Field, std::string> fd2str {
        {u, "u"}, {v, "v"}, {w, "w"}, {p, "p"}, {vertex, "vertex"}};

    // default values of type::str2bt
    std::map<std::string, BCType> str2bt {{"NOBC", NOBC},
        {"DIRICHLET", DIRICHLET}, {"NEUMANN", NEUMANN}, 
        {"CONVECTIVE", CONVECTIVE}, {"PERIODIC", PERIODIC}};
    // default values of type::bt2str
    std::map<BCType, std::string> bt2str {{NOBC, "NOBC"},
        {DIRICHLET, "DIRICHLET"}, {NEUMANN, "NEUMANN"}, 
        {CONVECTIVE, "CONVECTIVE"}, {PERIODIC, "PERIODIC"}};

    // default values of type::str2bl
    std::map<std::string, BCLoc> str2bl {
        {"left", XMINUS}, {"right", XPLUS}, {"bottom", YMINUS}, 
        {"top", YPLUS}, {"back", ZMINUS}, {"front", ZPLUS},
        {"xMinus", XMINUS}, {"xPlus", XPLUS}, {"yMinus", YMINUS}, 
        {"yPlus", YPLUS}, {"zMinus", ZMINUS}, {"zPlus", ZPLUS},
        {"XMINUS", XMINUS}, {"XPLUS", XPLUS}, {"YMINUS", YMINUS}, 
        {"YPLUS", YPLUS}, {"ZMINUS", ZMINUS}, {"ZPLUS", ZPLUS}};
    // default values of type::bl2str
    std::map<BCLoc, std::string> bl2str {
        {XMINUS, "xMinus"}, {XPLUS, "xPlus"}, {YMINUS, "yMinus"}, 
        {YPLUS, "yPlus"}, {ZMINUS, "zMinus"}, {ZPLUS, "zPlus"}};

} // end of namespace type
} // end of namespace petibm
