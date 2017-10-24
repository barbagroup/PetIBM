/*
 * io_singlebody.cpp
 * Copyright (C) 2017 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */


// STL
# include <fstream>
# include <sstream>

// PetIBM
# include <petibm/io.h>


namespace petibm
{
namespace io
{
    
PetscErrorCode readLagrangianPoints(const std::string &file, 
        PetscInt &nPts, type::RealVec2D &coords)
{
    PetscFunctionBeginUser;

    std::ifstream   inFile(file);
    std::string     line;
    PetscInt        c = 0;
    PetscInt        dim = 0;

    // check if we successfully open the file
    if (! inFile.good()) SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
            "Opening or reading body file %s failed!", file.c_str());


    // read the total number of points
    if (std::getline(inFile, line)) 
    {
        std::stringstream   sline(line);

        if (! (sline >> nPts))
            SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ, 
                    "Can't read the total number of points in file %s !\n",
                    file.c_str());

        if (sline.peek() != EOF) 
            SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ, 
                    "The first line in file %s contains more than one integer. "
                    "Please check the format.\n", file.c_str());
    }
    else
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
                "Error while reading the first line in file %s !\n", file.c_str());

    
    // get the dimension from the first coordinate set; initialize coords
    {
        std::getline(inFile, line);
        std::stringstream   sline(line);
        PetscReal temp;
        while (sline >> temp) dim += 1;
    
        if (dim == 0)
            SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
                    "Could not calculate the dimension from the first coordinate "
                    "set in the file %s!\n", file.c_str());

        // initialize the size of coordinate array
        coords = type::RealVec2D(nPts, type::RealVec1D(dim, 0.0));
        
        // read again to get first coordinate set
        sline = std::stringstream(line);
        for(PetscInt d=0; d<dim; ++d) sline >> coords[0][d];
        
        // increade c
        c += 1;
    }

    // loop through all points (from the second coordinate set)
    while(std::getline(inFile, line))
    {
        std::stringstream   sline(line);

        for(PetscInt d=0; d<dim; ++d)
            if (! (sline >> coords[c][d]))
                SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
                        "The number of doubles at line %d in file %s does not "
                        "match the dimension.\n", c+2, file.c_str());

        if (sline.peek() != EOF)
            SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
                        "The number of doubles at line %d in file %s does not "
                        "match the dimension.\n", c+2, file.c_str());

        c += 1;
    }

    // close the file
    inFile.close();

    // check the total number of points read
    if (c != nPts)
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ,
                "The total number of coordinates read in does not match the "
                "number specified at the first line in file %s !\n",
                file.c_str());

    PetscFunctionReturn(0);
}


PetscErrorCode print(const std::string &info)
{
    PetscFunctionBeginUser;
    
    PetscErrorCode  ierr;
    
    ierr = PetscSynchronizedPrintf(
            PETSC_COMM_WORLD, info.c_str()); CHKERRQ(ierr);
    
    ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT); CHKERRQ(ierr);
    
    ierr = PetscPrintf(PETSC_COMM_WORLD, "\n"); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

} // end of namespace io
} // petibm
