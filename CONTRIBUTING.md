# How to contribute to PetIBM

Welcome to the developer's guide of PetIBM!

## Adding new features and fixing bugs

All new features and bug fixes must go through a pull-request review procedure.
If you want to contribute to PetIBM, please fork the Barbagroup's [PetIBM](https://github.com/barbagroup/PetIBM) repository, make your changes on your fork, and then open a pull-request with your changes against the main PetIBM repository.

For new features and minor bugs (with small impact), the base branch of the pull-request should be the `develop` branch of the main repository.
(The `develop` branch will be merged into the `master` one once we are ready for a new release of PetIBM.)

For major bugs, the base branch should be the `master` branch of the main repository; it will be considered as a hotfix (bugfix) and a new version of PetIBM will be released as soon as possible by the maintainers with the micro number incremented.

New features should come with some kind of test or example to verify and/or validate the implementation.


## Reporting bugs and requesting new features

To report bugs, request new features, or simply ask questions, please open a GitHub issue on the Barbagroup's PetIBM repository.


## Writing documentation

New classes, methods, functions, and namespaces must be documented with Doxygen-style doctrings.

To locally generate and check the Doxygen documentation, use the command-line:

    > cd doc
    > doxygen Doxyfile

and open the file `doc/html/index.html` in your favorite browser.

You should also add code documentation whenever necessary; it will greatly help other developers to review your new features and bug fixes.

For new features, user's documentation must also be written.
For this purpose, we use Markdown files that are located in the `doc` folder of the root directory of PetIBM.

The Wiki pages and Doxygen API documentation are up-to-date with the latest release of PetIBM, which should also be the latest commit on the `master` branch.

Once a new release is drafted, we will merge the `doc` folder into the Wiki of PetIBM and the API documentation will be merged into the branch `gh-pages`.
