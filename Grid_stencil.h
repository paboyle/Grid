//////////////////////////////////////////////////////////////////////////////////////////
// Must not lose sight that goal is to be able to construct really efficient
// gather to a point stencil code. CSHIFT is not the best way, so probably need
// additional stencil support.
//
// Stencil based code could pre-exchange haloes and use a table lookup for neighbours
//
// Lattice <foo> could also allocate haloes which get used for stencil code.
//
// Grid could create a neighbour index table for a given stencil.
// Could also implement CovariantCshift.
//////////////////////////////////////////////////////////////////////////////////////////
