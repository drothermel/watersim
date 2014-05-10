#include "Water.h"
#include "MD.h"

int main(int argc, char *argv[])
{
	MD* moldyn = new MD(9.0, 4, 300, 100000, 0); //(boxL, nmols, T, steps, type)

	printf("Made it to the end!");
}