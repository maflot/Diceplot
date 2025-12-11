.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "NOTE: The 'diceplot' package is now legacy/maintenance-only.\n",
    "      For new projects, please use 'ggdiceplot': https://github.com/maflot/ggdiceplot \n\n",
    "      Or just install via ```install.packages(\"ggdiceplot\")```"
  )
}
