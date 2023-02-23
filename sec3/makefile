CPP = g++
CPPFLAGS = -O3 -Wall -Wextra -Wno-sign-compare -Wunused-function -fopenmp -std=c++14
OBJDIR = build/
SRCS = $(wildcard src/*.cpp)
OBJS = $(addprefix $(OBJDIR), $(SRCS:.cpp=.o))
DEPS = $(addprefix $(OBJDIR), $(SRCS:.cpp=.d))

$(OBJDIR)%.o: %.cpp
	$(info Generating $@)
	@mkdir -p $(dir $@)
	@$(CPP) -c $(CPPFLAGS) -o $@ $<

$(OBJDIR)%.d: %.cpp
	$(info Generating $@)
	@mkdir -p $(dir $@)
	@set -e; rm -f $@; \
	$(CPP) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,$(basename $(notdir $@)).o:,$(basename $@)\.o $(basename $@)\.d:,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

ifeq ($(filter main, $(MAKECMDGOALS)), main)
include $(DEPS)
endif

main: $(OBJS)
	$(CPP) $(CPPFLAGS) $^ -o $@
clean:
	-rm -r $(wildcard $(OBJDIR)*) main