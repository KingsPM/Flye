.PHONY: all clean debug profile

CXXFLAGS += -Wall -Wextra -pthread -std=c++11 -g -march=native
CFLAGS += -Wall -Wextra -march=native -g
LDFLAGS += -pthread -std=c++11 -rdynamic -march=native

ASSEMBLE_BIN := ${BIN_DIR}/flye-assemble
REPEAT_BIN := ${BIN_DIR}/flye-repeat
POLISH_BIN := ${BIN_DIR}/flye-polish

profile: CXXFLAGS += -pg
profile: LDFLAGS += -pg
profile: release

release: CXXFLAGS += -D_LOG -O3 -DNDEBUG
release: CFLAGS += -D_LOG -O2 -DNDEBUG
release: flye-repeat flye-polish flye-assemble

debug: CXXFLAGS += -Og -D_GLIBCXX_DEBUG -D_DEBUG -D_DEBUG_LOG -D_LOG
debug: CFLAGS += -Og -D_GLIBCXX_DEBUG -D_DEBUG -D_DEBUG_LOG -D_LOG
debug: flye-repeat flye-polish flye-assemble


#sequence module
sequence_obj := ${patsubst %.cpp,%.o,${wildcard sequence/*.cpp}} 

sequence/%.o: sequence/%.cpp sequence/%.c sequence/*.h common/*h
	${CXX} -c ${CXXFLAGS} $< -o $@

#ksw2 library
ksw2_obj := ${patsubst %.c,%.o,${wildcard ksw2/*.c}}

ksw2/%.o: ksw2/%.c ksw2/*h
	${CC} -c ${CFLAGS} $< -o $@

#flye-assemble module
assemble_obj := ${patsubst %.cpp,%.o,${wildcard assemble/*.cpp}}
flye-assemble: ${assemble_obj} ${sequence_obj} ${ksw2_obj}
	${CXX} ${assemble_obj} ${sequence_obj} ${ksw2_obj} -o ${ASSEMBLE_BIN} ${LDFLAGS}

assemble/%.o: assemble/%.cpp assemble/*.h sequence/*.h common/*h
	${CXX} -c ${CXXFLAGS} $< -o $@

#flye-repeat module
repeat_obj := ${patsubst %.cpp,%.o,${wildcard repeat_graph/*.cpp}}
flye-repeat: ${repeat_obj} ${sequence_obj} ${ksw2_obj}
	${CXX} ${repeat_obj} ${sequence_obj} ${ksw2_obj} -o ${REPEAT_BIN} ${LDFLAGS}

repeat_graph/%.o: repeat_graph/%.cpp repeat_graph/*.h sequence/*.h common/*.h
	${CXX} -c ${CXXFLAGS} $< -o $@

#flye-polish module
polish_obj := ${patsubst %.cpp,%.o,${wildcard polishing/*.cpp}}
flye-polish: ${polish_obj}
	${CXX} ${polish_obj} -o ${POLISH_BIN} ${LDFLAGS}

polishing/%.o: polishing/%.cpp polishing/*.h common/*h
	${CXX} -c ${CXXFLAGS} $< -o $@



clean:
	-rm ${repeat_obj}
	-rm ${ksw2_obj}
	-rm ${sequence_obj}
	-rm ${assemble_obj}
	-rm ${polish_obj}
	-rm ${REPEAT_BIN}
	-rm ${ASSEMBLE_BIN}
	-rm ${POLISH_BIN}
