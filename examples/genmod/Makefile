.SUFFIXES:
.SECONDARY:
.PHONY: all clean distclean help

mod1   := 2layer
mod2   := Viscoelastic_TTI_Thomsen

CC     := gcc
CFLAGS := -Wall -Wextra -O2 -g

help:
	@echo
	@echo -e "\033[4mPrograms to create simple exemplary models\033[0m"
	@echo
	@echo "Type '$(MAKE) all' to compile programs."
	@echo "Type '$(MAKE) clean' to clean directory."
	@echo

all: $(mod1) $(mod2)

$(mod1): $(mod1).c
$(mod2): $(mod2).c

%: %.c
	$(CC) $(CFLAGS) -o$@ $^ -lm

clean: 
	@rm -f $(mod1) $(mod2)
	@rm -f *~

distclean: clean
