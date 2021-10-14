.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "This is ragp >= 0.3.5 which adds several new features ",
    "and consequently breaking changes.\n",
    "Please read the NEWS: https://github.com/missuse/ragp/blob/master/NEWS.md.\n",
    "If you encounter any problems please report them: \n",
    "https://github.com/missuse/ragp/issues"
    )
}
