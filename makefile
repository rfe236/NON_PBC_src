.DEFAULT_GOAL := all

#include ${PETSC_DIR}/variables
#include ${PETSC_DIR}/rules

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

OBJS            = cubic_spline.o eam_setfl.o boundaryconditions.o input.o initialize.o minimizer.o neighbors.o output.o shells.o potential.o quadrature.o electrondensity.o calculateresults.o readresults.o main.o
#CFLAGS          = -I${BOOST_INC}
#CPPFLAGS        = -O2 -std=c++11
CLEANFILES      = DMD

all: $(OBJS) parser
	$(CLINKER) $(OBJS) -o DMD ${PETSC_TAO_LIB} -Lparser -lparser
#	${RM} ${OBJS}

main.o: main.cpp Vector3D.h input.h output.h minimizer.h cubic_spline.h eam_setfl.h
	$(CXX) -c main.cpp ${CXX_FLAGS} ${CXXFLAGS} ${CCPPFLAGS}

boundaryconditions.o: boundaryconditions.cpp Vector3D.h input.h
	$(CXX) -c boundaryconditions.cpp ${CXX_FLAGS} ${CXXFLAGS} ${CCPPFLAGS}

input.o: input.cpp input.h parser/Assigner.h parser/Dictionary.h
	$(CXX) -c input.cpp  ${CXX_FLAGS} ${CXXFLAGS} ${CCPPFLAGS}

initialize.o: initialize.cpp Vector3D.h input.h
	$(CXX) -c initialize.cpp  ${CXX_FLAGS} ${CXXFLAGS} ${CCPPFLAGS}

minimizer.o: minimizer.cpp minimizer.h Vector3D.h input.h potential.h KDTree.h
	${CXX} -c minimizer.cpp ${CXX_FLAGS} ${CXXFLAGS} ${CCPPFLAGS}

neighbors.o: neighbors.cpp Vector3D.h KDTree.h input.h
	$(CXX) -c neighbors.cpp ${CXX_FLAGS} ${CXXFLAGS} ${CCPPFLAGS}

output.o: output.cpp output.h input.h Vector3D.h 
	$(CXX) -c output.cpp ${CXX_FLAGS} ${CXXFLAGS} ${CCPPFLAGS}

shells.o: shells.cpp Vector3D.h
	$(CXX) -c shells.cpp ${CXX_FLAGS} ${CXXFLAGS} ${CCPPFLAGS}

potential.o: potential.cpp input.h potential.h cubic_spline.h eam_setfl.h
	$(CXX) -c potential.cpp ${CXX_FLAGS} ${CXXFLAGS} ${CCPPFLAGS}

quadrature.o: quadrature.cpp Vector3D.h
	$(CXX) -c quadrature.cpp ${CXX_FLAGS} ${CXXFLAGS} ${CCPPFLAGS}

electrondensity.o: electrondensity.cpp Vector3D.h potential.h
	$(CXX) -c electrondensity.cpp ${CXX_FLAGS} ${CXXFLAGS} ${CCPPFLAGS}

calculateresults.o: calculateresults.cpp Vector3D.h
	$(CXX) -c calculateresults.cpp ${CXX_FLAGS} ${CXXFLAGS} ${CCPPFLAGS}

readresults.o: readresults.cpp Vector3D.h input.h
	$(CXX) -c readresults.cpp ${CXX_FLAGS} ${CXXFLAGS} ${CCPPFLAGS}

eam_setfl.o: eam_setfl.cpp cubic_spline.h eam_setfl.h
	$(CXX) -c eam_setfl.cpp ${CXX_FLAGS} ${CXXFLAGS} ${CCPPFLAGS}

cubic_spline.o: cubic_spline.cpp cubic_spline.h
	$(CXX) -c cubic_spline.cpp ${CXX_FLAGS} ${CXXFLAGS} ${CCPPFLAGS}


parser::
##	(cd parser; make "CXX=$(CXX)")
	+make -C parser "CXX=$(CXX)"

allclean:
	make clean
	(cd parser; $(MAKE) clean)
