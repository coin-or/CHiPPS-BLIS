#==========================================================================#
# This file is part of the BiCePS Linear Integer Solver (BLIS).            #
#                                                                          #
# BLIS is distributed under the Eclipse Public License as part of the      #
# COIN-OR repository (http://www.coin-or.org).                             #
#                                                                          #
# Authors:                                                                 #
#                                                                          #
#          Yan Xu, Lehigh University                                       #
#          Ted Ralphs, Lehigh University                                   #
#                                                                          #
# Conceptual Design:                                                       #
#                                                                          #
#          Yan Xu, Lehigh University                                       #
#          Ted Ralphs, Lehigh University                                   #
#          Laszlo Ladanyi, IBM T.J. Watson Research Center                 #
#          Matthew Saltzman, Clemson University                            #
#                                                                          #
#                                                                          #
# Copyright (C) 2001-2023, Lehigh University, Yan Xu, and Ted Ralphs.      #
# All Rights Reserved.                                                     #
#==========================================================================#

## $Id$

include BuildTools/Makemain.inc

########################################################################
#                          Subdirectories                              #
########################################################################

SUBDIRS = src

pkgconfiglib_DATA = blis.pc

########################################################################
#             Additional files to be included in tarball               #
########################################################################

# Here we need include all files that are not mentioned in other Makefiles

EXTRA_DIST = \
	examples/BlisMain.cpp \
	examples/blis.par \
	examples/egout.mps \
	examples/flugpl.mps \
	examples/gen.mps \
	examples/knap1.mps \
	examples/knap2.mps \
	examples/knap3.mps \
	examples/VRP/AUTHORS \
	examples/VRP/INSTALL \
	examples/VRP/LICENSE \
	examples/VRP/Makefile.in \
	examples/VRP/README \
	examples/VRP/VrpCommonTypes.h \
	examples/VRP/VrpConstants.h \
	examples/VRP/VrpCutGenerator.cpp \
	examples/VRP/VrpCutGenerator.h \
	examples/VRP/VrpHeurTSP.cpp \
	examples/VRP/VrpHeurTSP.h \
	examples/VRP/VrpMain.cpp \
	examples/VRP/VrpModel.cpp \
	examples/VRP/VrpModel.h \
	examples/VRP/VrpNetwork.cpp \
	examples/VRP/VrpNetwork.h \
	examples/VRP/VrpParams.cpp \
	examples/VRP/VrpParams.h \
	examples/VRP/VrpSolution.h \
	examples/VRP/VrpVariable.h

########################################################################
#                           Extra Targets                              #
########################################################################

test: all
	cd test; $(MAKE) test

unitTest: test

clean-local: clean-doxygen-docs
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

install-exec-local: install-doc

uninstall-local: uninstall-doc uninstall-doxygen-docs

.PHONY: test unitTest doxydoc

