#Compilateur et option de compilation
CXX = g++
CXXFLAGS = -g -O3 -funroll-loops -march=native -Wall -Wextra -fsanitize=address -std=c++23

#Noms de fichier sources, objets et exécutable finale
SOURCES = main.cpp function.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = prog

#Règle principale pour compiler le projet
all: $(EXECUTABLE)

#Compilation de l'executable finale
$(EXECUTABLE) : $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(EXECUTABLE) $(OBJECTS)


#Compilation de chaque fichier .cpp en .o (objet)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

#Nettoyer les fichiers objets (.o) et l'executable

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

#Nettoyer tout et recompiler
rebuild: clean all