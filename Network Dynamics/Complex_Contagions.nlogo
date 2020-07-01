;;; QUESTION:  top of page 719 ...  "without replacement" ?? s

globals [
  k                                    ;; determines the average number of neighbors each node has
  this-cluster                         ;; a tool used by the identify-clusters procedure
  saturated?                           ;; When this is TRUE, no more agents can be activated
]

turtles-own [
  cluster                              ;; a tool used by the identify-clusters procedure
  distance-from-other-turtles          ;; list of distances of this node from other turtles, used by the identify-clusters procedure
  attempted?                           ;; keeps track of which turtles have state-checked since a change in the system
                                       ;; the attempted? variable is used to help the model be more efficient
]

links-own
[
  rewired?                             ;; keeps track of whether the link has been rewired or not, used by the rewire-all procedure
]





;;  The "setup" procedure sets up the model.
;;  First, it creates the network using "rewire-all."
;;  Then, it seeds the contagion

to setup
  clear-all
  ask patches [ set pcolor white ]

  set k 4
  rewire-all



  ask turtles [ set color black set attempted? false ]
  set saturated? false
  if count turtles > 0 [ seed-contagion ]
  reset-ticks
end



;;  The "construct-agent" procedure simply constructs an agent.
;;  In more complex models, a procedure like this
;;  can be helpful in generating agents with many different
;;  parameters or rules.

to construct-agent
  ;; turtles-own agent constructor
  set shape "circle"
  set label ""
end




;;  The "restart" procedure resets the contagion process without
;;  changing the network structure.

to restart
  ask turtles [ set color black set attempted? false ]
  set saturated? false
  if count turtles > 0 [ seed-contagion ]
  set-current-plot "success of contagion"
  clear-plot
  reset-ticks
end


;;  The "run-step" procedure runs a single step of the diffusion process.
;;  With each step, the model tries all the agents that are black (not red).
;;  When an agent is checked, it either updates (turns red) or remains the same.
;;  If the agent doesn't change, we set attempted? to TRUE for that turtle.

to run-step

  ;;  Check to see whether there are any turtles that might still become activated.
  ;;  If not, set "saturated?" to true and stop the function.
  ;;  If all the agents are red (none of them are black) then we know it's saturated
  ;;  We also know it's saturated if there are zero black agents who we haven't checked already.

  if count turtles with [ color = black ] = 0 [ set saturated? true stop ]
  if count turtles with [ color = black and attempted? = false ] = 0 [ set saturated? true stop ]

  if not saturated? [
    ask one-of turtles with [color = black ] [
      let count-triggered 0
      set count-triggered count link-neighbors with [color = red]

      ;;  If enough neighbors are active (red), an inactive (black) agent will be activated.
      ;;  Once a turtle turns from black to red, it's possible their neighbors may wish to change as well.
      ;;  Thus, whenever a turtle is activated, we reset attempted? to false for the other turtles.

      ifelse count-triggered >= threshold [
        set color red ask turtles with [ self != myself ] [ set attempted? false ]
      ]
      [
        set attempted? true
      ]

    ]
  ]
  tick
end




;;  This procedure "seeds" the contagion by randomly selecting a turtle
;;  and activating that turtle, as well as all of their neighbors (turning them red)

to seed-contagion
  ask one-of turtles [
    ask link-neighbors [ set color red ]
    set color red
  ]
end








;;  The "create-ringlat" procedure creates a new ring lattice.
;;  All this procedure does is create N turtles, set their color black,
;;  and then run the "wire-ringlat" procedure.

to create-ringlat
  crt N [ set color black construct-agent ]
  wire-ringlat
end






;;  The "wire-ringlat" procedure contains the machinery needed to wire a ring lattice.
;;  The parameter k determines the number of neighbors each node will have -- in this
;;  version, k is fixed.  However,  this parameter can be turned into a variable controlled
;;  from the interface.

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
    while [z <= floor (k / 2)]
    [
      ask turtle ni [ create-link-with turtle ((ni + z) mod count turtles) [ set rewired? false ] ]
      set z z + 1
    ]
    set ni ni + 1
  ]
end









;; the identify-clustesr and grow-clusters procedure ensure count the number of
;; components in the network.  This code is inspired by the NetLogo model "Dissemination of Culture"
;; by Iain Weaver.  You can find this model at: http://ccl.northwestern.edu/netlogo/models/community/Dissemination%20of%20Culture

to-report identify-clusters
  let max-cluster 0
  let num-clusters 0

  let seed one-of turtles
  ask turtles [ set cluster nobody ]
  while [seed != nobody] [
    ask seed [
      set cluster self
      set num-clusters num-clusters + 1
      set this-cluster 1
      grow-cluster
    ]
    if this-cluster > max-cluster [ set max-cluster this-cluster]
    set seed one-of turtles with [cluster = nobody]
  ]
  report list num-clusters max-cluster
end

to grow-cluster

    ask link-neighbors with [cluster = nobody] [
      if cluster = nobody [ set this-cluster this-cluster + 1 ]
      set cluster [cluster] of myself
      grow-cluster
    ]
end





;;  The "rewire-all" procedure generats a Small-World network according to the algorithm
;;  developed by Watts & Strogatz (1998).  This code is adapted from the Small Worlds NetLogo model
;;  developed by Uri Wilensky.  Original code Copyright 2005 Uri Wilensky.   See info tab for more details.

to rewire-all

  create-ringlat

  ;; set up a variable to see if the network is connected
  let success? false

  ;; if we end up with a disconnected network, we keep trying, because the APL distance
  ;; isn't meaningful for a disconnected network.
  let count-tries 0
  while [not success?] [
    ;; kill the old lattice, reset neighbors, and create new lattice
    ask links [ die ]
    wire-ringlat
;   set number-rewired 0

    ask links [

      ;; whether to rewire it or not?
      if (random-float 1) < p
      [
        ;; "a" remains the same
        let node1 end1
        ;; if "a" is not connected to everybody
        if [ count link-neighbors ] of end1 < (count turtles - 1)
        [
          ;; find a node distinct from node1 and not already a neighbor of node1
          let node2 one-of turtles with [ (self != node1) and (not link-neighbor? node1) ]
          ;; wire the new edge
          ask node1 [ create-link-with node2 [ set color cyan  set rewired? true ] ]

;          set number-rewired number-rewired + 1  ;; counter for number of rewirings
          set rewired? true
        ]
     ]
      ;; remove the old edge
      if (rewired?)
      [
        die
      ]
    ]

    set success? ( item 0 identify-clusters = 1 )
    set count-tries count-tries + 1
    if ( count-tries > 1000 ) [ set success? true print "couldn't make connected network!  try different parameters!" ]
  ]

end





;; This procedure rewires a single tie.
to rewire-one

  ;; make sure num-turtles is setup correctly else run setup first
  if count turtles != N [
    setup
  ]


  let potential-edges links with [ not rewired? ]
  ifelse any? potential-edges [
    ask one-of potential-edges [
      ;; "a" remains the same
      let node1 end1
      ;; if "a" is not connected to everybody
      if [ count link-neighbors ] of end1 < (count turtles - 1)
      [
        ;; find a node distinct from node1 and not already a neighbor of node1
        let node2 one-of turtles with [ (self != node1) and (not link-neighbor? node1) ]
        ;; wire the new edge
        ask node1 [ create-link-with node2 [ set color cyan  set rewired? true ] ]


        ;; remove the old edge
        die
      ]
    ]
  ]
  [ user-message "all edges have already been rewired once" ]
end





;;  This procedure counts the fraction (percent) of nodes who are activated (red)
to-report percent-saturated
  report ((count turtles with [color = red ] ) / (count turtles))
end
@#$#@#$#@
GRAPHICS-WINDOW
268
10
638
380
-1
-1
11.2
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

BUTTON
56
165
169
198
Setup Network
Setup
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
57
207
137
240
Go
ifelse not saturated? [ run-step ] [ stop ]
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
674
91
733
151
threshold
2
1
0
Number

PLOT
710
298
910
448
Success of Contagion
Time
% Saturated
0
10
0
1
true
false
"" ""
PENS
"% Saturation" 1 0 -16777216 true "" "plot percent-saturated"

BUTTON
57
255
196
288
Go Once
let totalblack count turtles with [ color = black ]\nwhile [ count turtles with [ color = black ] = totalblack and not saturated? ]\n[ run-step ]
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
840
165
941
202
Rewire One (SW)
rewire-one
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
12
10
252
142
Red nodes are \"active.\"   \n\nBlack nodes will activate if more than (threshold) of their neighbors are active.
20
0
1

TEXTBOX
838
35
1118
85
TOPOLOGY PARAMETERS
20
0
1

TEXTBOX
674
35
824
85
GAME PARAMETER
20
0
1

SLIDER
840
76
1030
109
p
p
0
1
0
0.01
1
NIL
HORIZONTAL

SLIDER
840
118
1012
151
N
n
10
100
56
2
1
NIL
HORIZONTAL

BUTTON
59
304
184
337
Restart Diffusion
restart
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1
@#$#@#$#@
## WHAT IS IT?

This model explores the spread of behavior through social networks. Each agent follows a simple threshold rule:  if enough neighbors adopt some behavior, they will adopt it as well.  When only one neighbor has to be activated, the  behavior is a simple contagion, and spread faster and more easily in small-world networks.  When the threshold is greater than 1/k, however, it is a complex contagion and even a small amount of randomness can disrupt the spread of adoption entirely.

The goal of this model is to understand how the structure of social networks can impact the spread of behavior.  A contagion can be any behavior that somebody might adopt, and the threshold captures the risk or uncertainty associated with that behavior:  the higher the uncertainty, the higher the threshold.

Under some network conditions, complex contagions spread very easily.  Under other conditions, however, some contagions will not spread at all!  Given two otherwise identical populations, we can see that the variation in the structure of interactions can have a big impact on social outcomes.

## HOW IT WORKS

The threshold t indicates the fraction of an agent's neighbors who must be active in order for a node to become active.  If an agent's has 5 neighbors and their t = 3/5, then they will be "active" if 3 or more neighbors are active, and "inactive" if 2 or fewer neighbors are active.

In this model, a contagion is "seeded" by activating a single node and all of its neighbors.  Then, the contagion spreads around the network!

Here's how the contagion spreads:  in each time-step, a randomly chosen non-active agent is selected.  That agent activates if their threshold is met. 

## HOW TO USE IT

Black nodes are inactive, and red nodes are active.

The THRESHOLD parameter determines the number of neighbors who must be active (red) in order for a node to activate (turn red).

The SETUP NETWORK button creates a "small world" network.  The TOPOLOGY parameters P and N determine the probability that a tie will be rewired and the number of nodes, respectively.  When P=0, a lattice network is created.  When a network is created, one node and all of their neighbors are randomly selected to be active (red).

The GO button runs the model until no more nodes can be activated.

The GO ONCE button randomly selects an eligible node, and actives them.

The RESTART DIFFUSION model resets the diffusion process by randomly selecting a new seed neighborhood without changing the network topology.

The REWIRE ONE button rewires a single tie.	


## THINGS TO NOTICE

Even slightly different network topoligies have substantively different collective behaviors.  For example, in a ring lattice with threshold 1/2 (0.5), even a single broken tie will stop a complex contagion:  try setting up a ring lattice and hitting the "rewire one" button once or a few times to see what happens.

How many ties can you rewire if THRESHOLD is less than 1/2?   How do these settings impact the time it takes for the contagion to spread through the whole network?


## THINGS TO TRY

Slowly raise THRESHOLD and see what happens.  What is the maximum threshold for a ring lattice?

Create a network with P=0 (a lattice network) and THRESHOLD=2.  Run the diffusion model a few times.  Now, rewire just a few times, and try the diffusion model a few times.  What happens?


## EXTENDING THE MODEL

What happens if threshold is not homogeneous?  

What happens if threshold is set so that nodes activated based on the fraction of active neighbors ("relative threshold") instead of a simple count ("absolute threshold")?


## RELATED MODELS

See the "Small Worlds With Diffusion" model that is available through the coursera course "Network Dynamics of Social Behavior" by Damon Centola.


## CREDITS AND REFERENCES


Centola, D., & Macy, M. (2007). Complex contagions and the weakness of long ties1. American Journal of Sociology, 113(3), 702-734.

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
