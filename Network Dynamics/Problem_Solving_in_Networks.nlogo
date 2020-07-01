;;;;
;;;;
;;;;    code created by
;;;;    joshua becker
;;;;    standing on the shoulders of giants
;;;;
;;;;



;;;;
;;;;  LEARNING POINTS
;;;;
;;;;  extension for this:
;;;;      when is the game over?
;;;;      adding the error factor
;;;;      allowing agents to change more than 1 bit at a time
;;;;      reducing the frequency with which agents copy
;;;;

extensions[array palette]
globals [
  nk-space         ;; an array of length B^K that determines the value of a given allele
  allele-map       ;; an array of length N that maps each digit-location to K other digit-locations
  infinity                             ;; a very large number.
                                         ;; used to denote distance between two turtles which
                                         ;; don't have a connected or unconnected path between them

  b

  global-max
  global-min

  plist            ;; used for
  numpeaks         ;; value distribution

]

turtles-own [
  solution                      ;; a turtle's solution to the complex problem.
  distance-from-other-turtles   ;; list of distances of this node from other turtles
                                ;; (just to make sure the small-world rewiring doesn't leave islands
]

links-own [
  rewired?
]

to setup
  clear-all
  ask patches [ set pcolor white ]
  set infinity 99999  ;; just an arbitrary choice for a large number
  set b 2

  ;; initialize distribution
  set numpeaks 8
  set plist []
  let i 0
  while [i < numpeaks]
  [
    set plist lput random-float 1 plist
    set i i + 1
  ]


  define-allele-values
  set-interdependencies

  spawn-turtles
  wire-ringlat
  reset-ticks
end

to go
  run-step
  let keepgoing false
  let compare-to [solution] of turtle 0
  ask turtles [
    ; If anybody's doesn't match, then we keep going.
    if (solution != compare-to) [ set keepgoing true ]
  ]
  if (not keepgoing) [ stop ]
end

;;;  run a step   ;;;

to run-step
  ask turtles [
    ;; step 1:  ask around.  do any of my neighbors have better solutions?
    let best-solution solution
    let better-one? false
    ask link-neighbors [
      if evaluate-fitness solution > evaluate-fitness best-solution
      [
        set best-solution solution
        set better-one? true
      ]
    ]

    ;; step 2:  if nobody had a better one, explore
    if not better-one?
    [
      let found-solution explore
      if evaluate-fitness found-solution > evaluate-fitness solution
      [
        set best-solution found-solution
        set better-one? true
      ]
    ]

    set solution best-solution

    let this-color precision evaluate-fitness solution 2
    ifelse enable-space-analysis
    [ set color palette:scale-scheme "Sequential" "Reds" 7 this-color 0 1 ]
    [ set label this-color ]
  ]
  tick
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;                                 ;;;
;;; the following set of functions  ;;;
;;; defines the agent behavior.     ;;;
;;;                                 ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


to spawn-turtles
  create-turtles Number-Of-Agents [
    ;; each turtle starts with a randomly generated solution
    set solution n-values n [ random b ]
    set label-color black
    set color palette:scale-scheme "Sequential" "Reds" 7 0 0 1
  ]
  set-default-shape turtles "circle"
  layout-circle (sort turtles) max-pxcor - 1
end


to-report explore
  ;; set the variable 'new answer' to be the turtle's solution variable
  ;; with a randomly chosen item (from 0 to n) replaced with a random value
  ;; from 0 to b that does not include the one already there.
  let replaced-item random n
  let new-answer replace-item replaced-item solution (item random (b - 1) list-not-x (item replaced-item solution) 0 b )
  report new-answer
end




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;                               ;;;
;;; the following set of methods  ;;;
;;; also known as 'functions'     ;;;
;;; also known as 'to-do-things'  ;;;
;;; defines the nk space.         ;;;
;;;                               ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;                                                                                           ;;
  ;; NK space is generated according to Appendix A in                                          ;;
  ;;      Lazer, D., & Friedman, A. (2007). The network structure of exploration and           ;;
  ;;      exploitation. Administrative Science Quarterly, 52(4), 667-694.                      ;;
  ;; There are two key components to this nk-space:                                            ;;
  ;;    a)  the value of a given "allele," which is here represented as a numeric string.      ;;
  ;;        - an agent has N alleles                                                           ;;
  ;;    b)  the pairs of digits that make up each allele.  these are randomly chosen           ;;
  ;;        - each allele has k digits                                                         ;;
  ;;                                                                                           ;;
  ;; The one departure from Lazer's implementation is that here we allow 'b' to be >2, meaning ;;
  ;; that our numeric strings CAN be bitstrings, as his were, but do not have to be.           ;;
  ;;                                                                                           ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




to define-allele-values
  ;; this function defines the value of each
  ;; allele of length k.  together with set-interdependencies
  ;; this defines the nk-space by creating an identifiable
  ;; value for any nk-string
  ;;
  ;; we're not actually directly defining the alleles:  rather,
  ;; we're defining them by their index in a list.  later,
  ;; a converter will take an allele (0 1 1) and convert it
  ;; to a decimal index.
  ;;
  ;; for a given index - allele - the value is randomly chosen
  ;; from a given distribution.

  set nk-space []
  let i 0
  repeat B ^ (K + 1) [
    set nk-space lput nk-distribution nk-space
    set i i + 1
  ]

end

to set-interdependencies
  ;; for each digit-location in the numeric string,
  ;; this function sets (k - 1) other digit-locations
  ;; to define alleles of size k.

  if k > (n - 1)
  [
    show "Error!  K cannot be greater than N-1.  Setting K to N-1."
    set k (n - 1)
  ]

  set allele-map []
  let i 0
  repeat n
  [
    set allele-map lput fput i n-of k (list-not-x i 0 n) allele-map
    set i i + 1
  ]
end



to-report nk-distribution
  ;; this function simply outputs a random number
  ;; from the appropriate distribution for the
  ;; desired nk-space.


  let p item (random numpeaks) plist

  let bin_ct 0
  repeat 1000 [
      if random-float 1 < p
      [
        set bin_ct bin_ct + 1
      ]
  ]
  report bin_ct
end


to-report evaluate-fitness [ test-solution ]
  let fitness 0
  let i 0
  while [i < n]
  [
    let this-allele []
    foreach item i allele-map
    [ [?1] ->
      set this-allele lput item ?1 test-solution this-allele
    ]
    set fitness fitness + item convert-allele this-allele nk-space
    set i i + 1
  ]
  set fitness (fitness / n)

  ; If max-value has been calculated,
  ; normalize the score.

  if (global-max != 0) [
    set fitness fitness / global-max
    set fitness fitness ^ 8
  ]

  report fitness
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;                               ;;;
;;; some generic, non-model       ;;;
;;; utilities.                    ;;;
;;;                               ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


to-report list-not-x [x x.min x.max]
  ;; a little utility to draw-without-replacement
  ;; from the range x.min to x.max excluding x.
  let i x.min
  let report-list []
  repeat (x.max - x.min)
  [
    if i != x
    [
      set report-list lput i report-list
    ]
    set i i + 1
  ]
  report report-list
end

to-report convert-allele [string]
  ;; this utility converts a base b string into
  ;; a decimal number - the index location of that
  ;; allele's value.
  ;; eg, when b = 2, this is a decimal-to-binary converter.
  let i 0
  let decimal-output 0
  repeat length string
  [
    set decimal-output decimal-output + ( ( b ^ (length string - i - 1) ) * item i string )
    set i i + 1
  ]
  report decimal-output
end

to-report convert-to-allele [ decimal ]
  ;; reverse process of above - so if b=2, this
  ;; is a decimal-to-binary converter.


  let output-string []
  repeat n [
    set output-string fput (decimal mod b) output-string
    set decimal floor (decimal / b)
  ]

  report output-string
end











;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;
;;;
;;; much of this network construction code was
;;; copy-pasted from the uri wilensky small worlds model
;;;
;;; some of it was copy-pasted from previous work of mine
;;;
;;;
;;; remember:
;;;
;;; if you're not copy-pasting,
;;; you're not coding!
;;;
;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;



to wire-ringlat
  layout-circle (sort turtles) max-pxcor - 1

  layout-circle (sort turtles with [ who mod 2 = 0] ) max-pxcor - 4


  ;; iterate over the turtles
  let ni 0
  while [ni < count turtles]
  [
    ;; make edges with the next two neighbors
    ;; this makes a lattice with average degree of 4
    let z 1
    while [z <= floor (degree / 2)]
    [
      ask turtle ni [ create-link-with turtle ((ni + z) mod count turtles) ]
      set z z + 1
    ]
    set ni ni + 1
  ]
end



;; connects the two turtles
to make-edge [node1 node2]
  ask node1 [ create-link-with node2  [
    set rewired? false
  ] ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;                            ;;;
;;;   game state computations  ;;;
;;;                            ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report number-of-unique-solutions
  let solution-set []
  ask turtles
  [
    let already-there? false
    foreach solution-set
    [ [?1] ->
      if solution = ?1
      [
        set already-there? true
      ]
    ]
    if not already-there?
    [
      set solution-set lput solution solution-set
    ]
  ]

  report length solution-set
end


to-report average-score
  let avg 0
  ask turtles [
    set avg avg + evaluate-fitness solution
  ]
  report ( avg / count turtles )
end

to determine-peaks
  if enable-space-analysis [
    ;; WARNING!!!
    ;; this function can easily take a long time
    ;; if n or b are too high...

    ;; determine decimal value of highest numerical-string
    let mylist []
    repeat n
    [
      set mylist lput (b - 1) mylist
    ]

    let best-answer 0
    let worst-answer 999999999
    let i 0
    while [i <= convert-allele mylist]
    [
      let this-value evaluate-fitness (convert-to-allele i)
      if this-value > best-answer
      [
        set best-answer this-value
      ]
      if this-value < worst-answer
      [
        set worst-answer this-value
      ]
      set i i + 1
    ]

    set global-max best-answer
    set global-min worst-answer
  ]
end

to-report global-peak

  ;; report global max without recalculating
  report global-max
end

to-report global-minimum

  ;; report global min without recalculating

  report global-min
end


@#$#@#$#@
GRAPHICS-WINDOW
162
12
591
441
-1
-1
13
1
10
1
1
1
0
0
0
1
-16
16
-16
16
0
0
1
ticks
30

INPUTBOX
624
72
674
132
n
10
1
0
Number

INPUTBOX
679
72
729
132
k
5
1
0
Number

BUTTON
27
24
109
57
Setup
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
631
233
737
293
number-of-agents
40
1
0
Number

TEXTBOX
626
23
776
73
NK Space Variables
20
0
1

TEXTBOX
630
203
780
228
Network Setup
20
0
1

BUTTON
27
61
110
94
One Step
run-step
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
26
98
112
139
Go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
854
301
1166
459
Number of Unique Solutions
Time
Unique Solutions
0
20
0
20
true
false
"" ""
PENS
"default" 1 0 -16777216 true "" "plot number-of-unique-solutions"

PLOT
853
40
1184
190
Average Score
Time
Score
0
20
0
1000
true
false
"determine-peaks\nif enable-space-analysis [\n   set-plot-y-range 0 1\n]" ""
PENS
"default" 1 0 -16777216 true "" "plot average-score"

TEXTBOX
625
141
782
174
The nk-space is randomly generated every time you click 'setup'
9
0
1

SWITCH
854
203
1032
236
enable-space-analysis
enable-space-analysis
0
1
-1000

TEXTBOX
1039
203
1189
236
Warning!  NK Space analysis algorithm is very slow for large N and K.
9
0
1

SLIDER
629
308
801
341
degree
degree
2
50
4
2
1
NIL
HORIZONTAL
@#$#@#$#@
## WHAT IS IT?

This model is designed to demonstrate the effect of network topology on complex problem solving.  It is based on:

Lazer, D., & Friedman, A. (2007). The network structure of exploration and exploitation. Administrative Science Quarterly, 52(4), 667-694.

## HOW IT WORKS

Each agent starts with a randomly chosen "solution," represented by a bit string where each slot has B possible values.

At each time step, an agent first checks to see if any neighbors have a better solution.  If they do, the agent will "exploit" that social information by adopting the best solution from among their neighbors.

If neighbors do not offer a better solution, an agent will "explore" by randomly changing one bit of their bit-string solution; if that new solution is better, they keep that new solution.  The value of a solution is determined by a randomly generated NK space.


The exact details of the "NK Space" are outside of the scope of this description.  Suffice to say, an "NK Space" is a way of mapping solutions to complex problems on a fitness landscale offering "tunable ruggedness."  A "fitness landscape" is a landscape of all the possible solutions to a problem, where "fitness" refers to how good a solution is.

By increasing "K" you increase the interdependence between parts of the solution -- so that with higher K, changing one choice impacts the fitness of other choices.  (Think:  adding extra thrusters to a plane also increases weight, so the two decisions are interdependent.) 


## HOW TO USE IT

Set the network topology using the "topology" parameters.

Set the NK space parameters by setting N and K.  N is the length of the bit-string (the number of factors that must be considered in the solution) and K is the extent of interdependency, or complexity of the problem.

## THINGS TO NOTICE

For very simple problems (N=10, k=0) Who will find a better solution in the long run - agents in a fully connected network or agents in a line or lattice?  

What about for very complex problems (N=10, k=5)?

What prevents an agent or network from finding the optimal solution?

## EXTENDING THE MODEL

In this model, agents never make mistakes.  What happens if you add "noise," or random variations, to the model?

Agents here are limited ot changing one feature of their bitstring at a time.  What happens if they can change more?

## NETLOGO FEATURES

This model uses the 'array' extension to generate the NK space.

## CREDITS AND REFERENCES

Lazer, D., & Friedman, A. (2007). The network structure of exploration and exploitation. Administrative Science Quarterly, 52(4), 667-694.

## COPYRIGHT AND LICENSE

Copyright 2017 Joshua Becker

![CC BY-NC-SA 3.0](https://licensebuttons.net/l/by-nc-sa/3.0/88x31.png)

This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 License.  To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 559 Nathan Abbott Way, Stanford, California 94305, USA.

<!-- 2017 -->
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0
-0.2 0 0 1
0 1 1 0
0.2 0 0 1
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@

@#$#@#$#@
