#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

extern char const * program_name;

void set_signal_handler();
void posix_print_stack_trace();
static uint8_t alternate_stack[SIGSTKSZ];
int  divide_by_zero();
void cause_segfault();
void stack_overflow();
void infinite_loop();
void illegal_instruction();
void cause_calamity();
