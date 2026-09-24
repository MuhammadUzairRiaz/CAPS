// Process-wide settings the Studio changes (Settings board): the colour palette and the worker-thread cap.
#pragma once

namespace caps {

enum class Palette { Caps = 0, OkabeIto = 1, Monochrome = 2 };
void set_palette(Palette p);
Palette palette();
// Display colour of an element and of the k-th molecule / chain / type in the current palette (0xRRGGBB).
unsigned element_colour(int z);
unsigned category_colour(int k);

// Worker threads for the parallel loops (pair terms, Pack, analysis): 0 = one per hardware thread, at most 16.
void set_max_threads(int n);
int max_threads();

}  // namespace caps
