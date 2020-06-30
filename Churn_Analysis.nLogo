extensions [ nw ]

globals [
  number-quality-issues
]

turtles-own [

  churn? ;; identify if the customer has churned
  contract? ;; identify customers who have contracts with a given company
  remaining-ctr-time ;; identify the time to the end of contracts with a given company
  satisfaction-level
  quality-sensitiveness ;; scale from zero to one defining how much a change in price affects the customer perception of value
  price-sensitiveness ;; scale from zero to one defining how much a change in quality affects the customer perception of value
  network-sensitiveness ;; scale from zero to one defining how much a change in the network opinion affects the customer perception of value

  price-assessment ;; identify how the customer asses the price of the service
  quality-assessment ;; identify how the customer asses the quality of the service
  network-assessment ;; identify how the network of the customer assess the service

]

to calculate-satisfaction

  ;; in this model the network assessment on the quality is based on the direct links of an agent

  if any? link-neighbors [
    set network-assessment sum [satisfaction-level] of link-neighbors / count link-neighbors
  ]

  ;; the individual satisfaction is the average weight of quality, price and network influence

  let weighted-satisfaction quality-sensitiveness * quality-assessment + price-sensitiveness * price-assessment + network-sensitiveness * network-assessment
  let overall-weight quality-sensitiveness + price-sensitiveness + network-sensitiveness
  set satisfaction-level weighted-satisfaction / overall-weight

end

to setup-network-with-randomness

  ;; for this setup the network is initialized adding randomness to the parameters defined on the setup, in such a way that agents can have different parameters
  clear-all

  nw:generate-random turtles links number-of-agents connection-probability

  set-default-shape turtles "circle"

  let root-agent max-one-of turtles [ count my-links ]
  layout-radial turtles links root-agent

  set number-quality-issues 0
  ask turtles [

    set churn? False
    ifelse random-float 1.0 > contract-percentage
    [
      set contract? True
      set remaining-ctr-time random 12
      if remaining-ctr-time = 0 [
        set remaining-ctr-time 1
      ]
      set color green
    ]
    [
      set contract? False
      set remaining-ctr-time 0
    ]

    set color blue

    set quality-sensitiveness (quality-importance + random-float (.1 * quality-importance))
    set price-sensitiveness (price-importance + random-float (.1 * price-importance))
    set network-sensitiveness (network-importance + random-float (.1 * network-importance))

    if quality-sensitiveness > 1 [set quality-sensitiveness 1]
    if price-sensitiveness > 1 [set price-sensitiveness 1]
    if network-sensitiveness > 1 [set network-sensitiveness 1]

    ;; for increasing the randomness we consider that some customer will have opinions below the average and others above the average
    ifelse random-float 1.0 < .5
    [
      set price-assessment price-satisfaction-level - random-float(.1 * price-satisfaction-level)
      set quality-assessment service-satisfaction-level - random-float(.1 * service-satisfaction-level)
    ]
    [
      set price-assessment price-satisfaction-level + random-float(.1 * price-satisfaction-level)
      set quality-assessment service-satisfaction-level + random-float(.1 * service-satisfaction-level)
    ]

    if price-assessment < 0 [set price-assessment 0]
    if price-assessment > 1 [set price-assessment 1]

    if quality-assessment < 0 [set quality-assessment 0]
    if quality-assessment > 1 [set quality-assessment 1]

    set network-assessment 1 ;; for simplicity we are going to setup as initial value for the network assessment the value 1

    calculate-satisfaction

    set label int (satisfaction-level * 100)
  ]

  reset-ticks
end


to setup-network

  ;; for this setup the network is initialized with all agents have the same parameters defined on the setup
  clear-all

  nw:generate-random turtles links number-of-agents connection-probability

  set-default-shape turtles "circle"

  let root-agent max-one-of turtles [ count my-links ]
  layout-radial turtles links root-agent

  set number-quality-issues 0

  ask turtles [

    set churn? False
    ifelse random-float 1.0 > contract-percentage
    [
      set contract? True
      set remaining-ctr-time random 12
      if remaining-ctr-time = 0 [
        set remaining-ctr-time 1
      ]
      set color green
    ]
    [
      set contract? False
      set remaining-ctr-time 0
    ]

    set color blue

    set quality-sensitiveness quality-importance
    set price-sensitiveness price-importance
    set network-sensitiveness network-importance

    ;; for increasing the randomness we consider that some customer will have opinions below the average and others above the average
    ifelse random-float 1.0 < .5
    [
      set price-assessment price-satisfaction-level
      set quality-assessment service-satisfaction-level
    ]
    [
      set price-assessment price-satisfaction-level
      set quality-assessment service-satisfaction-level
    ]


    set network-assessment 1 ;; for simplicity we are going to setup as initial value for the network assessment the value 1

    calculate-satisfaction

    set label int (satisfaction-level * 100)
  ]

  reset-ticks
end



to go

  if all? turtles [ churn? = true ] [ stop ]

  ask turtles [

    calculate-satisfaction

    set label int (satisfaction-level * 100)

    ifelse contract? = True
    [
      set remaining-ctr-time (remaining-ctr-time - 1) ;; at each tick the contract duration will be reduced up to the limit of 0
      if remaining-ctr-time = 0 [
        ifelse satisfaction-level >=  threshold-contract-reneweal
        [
          set remaining-ctr-time 12
        ]
        [set contract? False]
      ]
    ]
    [
      if satisfaction-level < .6 [
        set churn? True
        set color red]
    ]
  ]

  tick
end

to quality-issue

  ;; quality issue will be applied to the percentage of customer affected defined as the parameter of the model

  ask turtles [
    if color != red [
      ifelse random-float 1.0 <= percentage-affected [
        if quality-assessment > 0 [
          set quality-assessment quality-assessment - .1
          calculate-satisfaction
          set color yellow
          set number-quality-issues number-quality-issues + 1
          set label int (satisfaction-level * 100)
        ]
      ]
      [set color blue]
    ]
  ]

end
@#$#@#$#@
GRAPHICS-WINDOW
229
49
979
800
-1
-1
21.2
1
12
1
1
1
0
0
0
1
-17
17
-17
17
1
1
1
ticks
30.0

BUTTON
19
552
198
585
Setup network
setup-network
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
490
10
846
48
Churn Analysis on Social Networks
15
0.0
0

PLOT
992
51
1494
380
Satisfaction level of Network
NIL
NIL
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -955883 true "" "plot sum [satisfaction-level] of turtles / count turtles"

PLOT
993
396
1495
731
Churn
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -2674135 true "" "plot count turtles with [churn? = true]"

BUTTON
19
636
196
669
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

TEXTBOX
20
680
170
708
Interactions between company and agents
11
75.0
1

SLIDER
18
87
196
120
number-of-agents
number-of-agents
1
100
50.0
1
1
NIL
HORIZONTAL

SLIDER
18
130
197
163
connection-probability
connection-probability
0
1
0.1
.01
1
NIL
HORIZONTAL

SLIDER
18
244
197
277
contract-percentage
contract-percentage
0
1
0.45
.01
1
NIL
HORIZONTAL

MONITOR
999
751
1169
796
agents with active contracts
count turtles with [contract? = True]
17
1
11

BUTTON
16
760
188
793
Report problem of quality
quality-issue
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
16
718
188
751
percentage-affected
percentage-affected
0
1
0.34
.01
1
NIL
HORIZONTAL

SLIDER
18
288
197
321
price-satisfaction-level
price-satisfaction-level
0
1
1.0
.01
1
NIL
HORIZONTAL

SLIDER
18
172
197
205
threshold-contract-reneweal
threshold-contract-reneweal
.5
1
0.7
.01
1
NIL
HORIZONTAL

SLIDER
17
332
196
365
service-satisfaction-level
service-satisfaction-level
0
1
0.54
.01
1
NIL
HORIZONTAL

TEXTBOX
23
44
191
62
Setup Parameters
12
0.0
1

SLIDER
20
419
198
452
price-importance
price-importance
0
1
0.73
.01
1
NIL
HORIZONTAL

SLIDER
20
461
198
494
quality-importance
quality-importance
0
1
0.65
0.01
1
NIL
HORIZONTAL

TEXTBOX
24
70
174
88
Network
11
5.0
1

TEXTBOX
21
221
171
239
Company statistical data
11
5.0
1

TEXTBOX
21
383
171
411
Customer decision variables importance
11
5.0
1

SLIDER
20
504
198
537
network-importance
network-importance
0
1
0.64
.01
1
NIL
HORIZONTAL

MONITOR
1181
750
1335
795
satisfaction level
(sum [satisfaction-level] of turtles / count turtles) * 100
17
1
11

BUTTON
17
590
197
623
Setup network with randomness
setup-network-with-randomness
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
1350
749
1497
794
number of quality issues
number-quality-issues
17
1
11

@#$#@#$#@
## WHAT IS IT?

This model intends to explore how quality issues for some individuals on a network can increase the risk of churn thanks to influential factors.

## HOW IT WORKS

When the model is set up it creates the number of individuals defined by the user in the *number-of-agents* and defines the links according the parameter *connection-probability*.

### Agents

The agents represent the individuals of the social network being created based on the parameter *number-of-agents*. The number of connections is defined based on the *connection-probability*. Everyone has its own characteristics that are represented by:

- **churn?** identify if the customer has churned. It is initialized as False for the whole network.

- **contract?** identify customers who have contracts with a given company. Customers with contracts usually tend to go up to the end of it to avoid fines, so we consider that customer will not churn before the termination of a contract. The number of agents defined with contracts is calculated based on the *contract-percentage*.

	The parameter *threshold-contract-reneweal* defines what should be the satisfaction level to a customer renewal a contract.

- **remaining-ctr-time** identify the time to the end of contracts with a given company. For agents with *contract?* equals True, it is initialized based on random number on the interval [0, 12].

- **quality-sensitiveness** scale from zero to one defining how much a change in price affects the customer perception of value. It is initialized based on random number on the interval [0, 1].

- **price-sensitiveness** scale from zero to one defining how much a change in quality affects the customer perception of value. It is initialized based on random number on the interval [0, 1].

- **network-sensitiveness** ;; scale from zero to one defining how much a change in the network opinion affects the customer perception of value.  It is initialized based on random number on the interval [0, 1].


- **price-assessment** identify how the customer asses the price of the service. 

- **quality-assessment** identify how the customer asses the quality of the service

- **network-assessment** identify how the network of the customer assess the service

- **satisfaction-level** it is calculated based on an equation that considers the average weight for satisfaction based on three elements price, quality and network influence.

### System Dynamics

This model works with the inertia concept, in a way that we need to have an external factor to provoke a change in the satisfaction level of an individual. These perturbations to the system are introduced when the user reports a problem of quality that affects agents according *percentage-affected*.

	ask turtles [
    if color != red [
      ifelse random-float 1.0 < percentage-affected [
        if quality-assessment > 0 [
          set quality-assessment quality-assessment - .1
          calculate-satisfaction
          set color yellow
          set number-quality-issues number-quality-issues + 1
          set label int (satisfaction-level * 100)
        ]
      ]
      [set color blue]
    ]  ]

When the model is running at each tick the agent’s satisfaction is updated and if the satisfaction level is lower than .6 and the customer do not have a contract then it will churn. At each iteration the *remaining-ctr-time* is decreased by one unit of time.

The formula used to calculate the satisfaction of an individual is a weighted average of price, quality and network influence perception. The level of importance that is used to weight is defined by the parameters that can be modified by the user: *price-importance*, *quality-importance* and *network-importance*.

	to calculate-satisfaction
    
  	if any? link-neighbors [
    	set network-assessment sum [satisfaction-level] of link-neighbors / count 	link-neighbors]
  	let weighted-satisfaction quality-sensitiveness * quality-assessment + price-sensitiveness * price-assessment + network-sensitiveness * network-assessment
	  let overall-weight quality-sensitiveness + price-sensitiveness + network-sensitiveness
 	 set satisfaction-level weighted-satisfaction / overall-weight
  
	end


### Assumptions

- the company will not invest in any marketing campaigns to raise the value perception of the agents.
- the model stops if all the customers have churned.


## HOW TO USE IT

1) Define the parameters you want to simulate moving the sliders to the values you identify as appropriate to your initial analysis. Then you need to set up the network and to do so you have two alternatives:

- SETUP NETWORK: in this case all agents will be initialized with the same properties that define the customer satisfaction based on the parameters of price, quality and network influence defined on the sliders.

- SETUP NETWORK WITH RANDOMNESS: here we introduce some randomness on the parameters, taking a standard deviation of 10% of the original parameter. This will allow that some agents be under and others above the average.

2) Press the button GO to start running the model and reach the stabilization of the network satisfaction level.

3) Introduce some perturbations using the button REPORT PROBLEM OF QUALITY, it will reduce the level of satisfaction of individuals and impact also on its neighbors.

## THINGS TO NOTICE

### Network Assessment

At the beginning we need some ticks to stabilize the overall satisfaction level , this happens due to the fact that we initialize all agents with the property *network-assessment* value set to 1, but at each tick we update the calculation based on the neighbors, as each individual has different parameters that compose the satisfaction we need some ticks to reach the actual initial state.

### Inertia

If no problem is reported there will be no events of churn. It means that our model is not affected by any other factor outside the context of customer and company.

### Churn

Check how many issues need to be reported to have a complete churn on the network for different configurations of the network.


## THINGS TO TRY

Try increasing the network-importance and check the influence on the churn when compared to smaller values for this parameter.

In the company statistical data simulate different configurations of customer level of satisfaction with price and quality.

## EXTENDING THE MODEL

Try modelling a scenario with competitors where we also simulate what would be the new distribution of the network considering that customers can move from one company to another when they are not satisfied.

Analyze a network created with real world data collected from companies using the idea of the model 'Network Import Example'.

Introduce some marketing initiatives that can improve the assessment of the individual and of the whole network.

## NETLOGO FEATURES

The primitive generate-random from NetLogo nw extension was used to easily create a random that creates links based on a predefined connection probability.

We used the agentset feature to identify the neighbors of a given individual and collecting the opinion of them on the service.

## RELATED MODELS

The model 'Simple Viral Marketing' allows you to explore the relationship between different centrality measures, and different network types, and see if the interactions between them make for faster or slower spread of a product.

## CREDITS AND REFERENCES

This model was inspired by the concepts of CUSTOMER-INFLUENCE FACTOR MODELING, in which the idea of integrating traditional machine learning techniques to agent based model in order to get better insights for business purposes, presented by:

Reis Pinheiro, Carlos Andre. Social Network Analysis in Telecommunications. Hoboken, New Jersey, Wiley, May 2011.

The network generation code was inspired by the library model 'NW General Examples', that demonstrates various features of the Netlogo NW extension.

## COPYRIGHT AND LICENSE

Copyright 2020 Erika Fonseca.

![CC BY-NC-SA 4.0](https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png)

This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.  To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to Creative Commons, 559 Nathan Abbott Way, Stanford, California 94305, USA.

Commercial licenses are also available. To inquire about commercial licenses, please contact Erika Fonseca at erika.luiz-fonseca@edu.dsti.institute.
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

line-half
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

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.1
@#$#@#$#@
import-network
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
