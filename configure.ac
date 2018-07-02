AC_PREREQ([2.69])
AC_INIT(intercrop, 0.1, expiredpopsicle@gmail.com)
AM_INIT_AUTOMAKE([subdir-objects])

AC_PROG_CXX
AC_PROG_RANLIB

PKG_CHECK_MODULES(libjpeg, libjpeg, [], [])
PKG_CHECK_MODULES(libpng,  libpng,  [], [])

liblilyengine_LIBS=$(Lily-Engine-Utils/lilyutils-config --libs)
liblilyengine_CFLAGS=$(Lily-Engine-Utils/lilyutils-config --cxxflags)
AC_SUBST(liblilyengine_LIBS)
AC_SUBST(liblilyengine_CFLAGS)

AC_OUTPUT(Makefile src/Makefile)
