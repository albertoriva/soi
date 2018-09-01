# soi

SOI (Signficance Of Intersections) is a command-line program to quickly evaluate whether the intersections among sets (up to four) are significantly high or low. It works by simulating a large number of random sets with the same number of elements as the original sets and measuring their interesections. The observed intersections are then compared against the distributions of random ones obtained through simulation in order to compute an empirical P-value.

## Usage

The syntax for this command is:

```
soi [-h] [-t] [-i I] [-n N] [-2] specs...
```

`specs` consists of entries of the form `S=n`, where S describes a set or an intersection, and n is the
number of elements it contains. S can be one of: N, A, B, C, D, AB, AC, AD, BC, BD, CD, ABC, ABD, ACD,
BCD, ABCD. N indicates the total number of items in the domain (can also be set with the `-n` option), 
single letters indicate the size of sets, multiple letters indicate intersections.

In 2-way mode, activated with the `-2` option, only two sets (A, B) and their intersection (AB) are used. 
The output consists of the expected and observed values for the intersection, for A only (A0) and for B 
only (B0) with the respective P-values. P(+) indicates that the observed value is higher than the expected one, while
P(-) indicates that it is lower.

If the `-r` option is supplied, the program will write an R script with the specified filename that, when executed,
will draw Venn diagrams of the sets and their intersections to file `venn.png`. The name of this file can
be changed with the `-i` option.

## Options

Option|Description
------|-----------
 -h         | Print this help message.
 -t         | Output results in tab-delimited format.
 -n N       | Set the total number of items to N (default: 1000).
 -i I       | Set the number of iterations to I (default: 100000).
 -r rfile   | Write an R script to plot the specified sets to `rfile`.
 -i imgfile | The R script will save image to png file `imgfile` (default: venn.png).
 -2         | Enable 2-way mode. Only uses the following specs: A, B, AB.

## Example

Imagine we conducted a survey of 1,000 individuals, and we found that 200 of them like
apples (A), 100 like bananas (B), and 50 like cherries (C). We also found that 22 or them 
like both apples and bananas (AB), 20 like both apples and cherries (AC), and 40 like 
bananas and cherries (BC). Finally, 5 individuals like all three fruits. We can ask the
question: are these intersections larger than we would expect by chance? E.g., is there
reason to believe that people who like apples also like bananas?

This can be answered with soi as follows:

```
$ soi N=1000 A=200 B=100 C=50 AB=22 AC=20 BC=40 ABC=5
# N = 1000
AB: Expected=20, Observed=22, P=0.249757
AC: Expected=10, Observed=20, P=0.000170
BC: Expected=5, Observed=40, P=0.000010
ABC: Expected=1, Observed=5, P=0.002950
```

Assuming a standard significance level of 1%, This result tells us that the intersection
between apple-likers and banana-likers is *not* significant. On the other hand, the
intersection between A and C and the one between B and C are highly significant. Finally,
the intersection of all three groups is also significant, although at a lower level.
