\ -- Variable, Constant, Array -----------------------------------------

: variable create 0 , does> ;
: constant create , does> @ ;
: ?        @ . ;
: c?       c@ . ;
: +!       dup @ rot + swap ! ;

\ -- Boolean constants -------------------------------------------------
-1 constant true
0  constant false

\ -- Stack helpers -----------------------------------------------------
: 2dup   over over ;
: 2drop  drop drop ;
: nip    swap drop ;
: tuck   swap over ;
: ?dup   dup if dup then ;

\ -- Arithmetic helpers ------------------------------------------------
: 1+     1 + ;
: 1-     1 - ;
: 2+     2 + ;
: 2-     2 - ;
: 2*     2 * ;
: 2/     2 / ;
: negate 0 swap - ;
: abs    dup 0 < if negate then ;
: min    2dup < if drop else swap drop then ;
: max    2dup > if drop else swap drop then ;

\ -- Logic -------------------------------------------------------------
\ Note: and/or/xor/invert are bitwise primitives.
\ not, <>, <=, >= are defined in terms of comparisons.
: not  0= ;
: <>   = not ;
: <=   > not ;
: >=   < not ;

\ -- Output helpers ----------------------------------------------------
: space    32 emit ;
: spaces   dup 0 > if 0 do space loop else drop then ;
: .cr      . cr ;

