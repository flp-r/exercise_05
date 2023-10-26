# Hirschberg's algorithm
Hirschberg's algorithm is a dynamic programming algorithm that finds the optimal sequence alignment between two strings
**X** with length *m* and **Y** with length *n*.
In comparison with Needleman-Wunsch algorithm, Hirschberg's algorithm is more space efficient:

| Algorithm | Time complexity | Space complexity |
|---|-----------------|------------------|
| Needleman-Wunsch | *O*(*m* \* *n*) | *O*(*m* \* *n*)  |
| Hirschberg | *O*(*m* \* *n*) | *O*(*m* + *n*)   |

### Task
Implement Hirschberg's algorithm - you can use a template `hirschberg_template.R`.


More information can be also find:
* in the presentation by Kevin Wayne [here](https://www.cs.princeton.edu/~wayne/kleinberg-tardos/pdf/06DynamicProgrammingII.pdf),
* in the presentation by Carl Kingsford [here](http://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/linspace.pdf),
* or on Wikipedia [page](https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm).

<details>
<summary>Spoilers! Help</summary>

#### Pseudo code of Hirschberg's algorithm
```
Hirschberg(X, Y, align, match, mismatch, gap)
1   Z ← the first row of alignment
2   W ← the second row of alignment
3   if length(X) = 0
4     for i ← 1 to length(Y)
5       Z ← Z + '-'
6       W ← W + Y[i]
7     align ← merge alignments Z and W
8   else if length(Y) = 0
9     for i ← 1 to length(X)
10      Z ← Z + X[i]
11      W ← W + '-'
12    align ← merge alignments Z and W
13  else if length(X) = 1 and length(Y) = 1
14    Z ← Z + X[1]
15    W ← W + Y[1]
16    align ← merge alignments Z and W
17  else
18    xlen ← length(X)
19    xmid ← xlen / 2
20    ylen ← length(Y)
21    ScoreL ← NWScore(X[1:xmid], Y, match, mismatch, gap)
22    ScoreR ← NWScore(reverse(X[xmid+1:xlen]), reverse(Y))
23    ymid ← arg max (ScoreL + reverse(ScoreR)) - 1
24    align ← Hirschberg(X[1:xmid], Y[1:ymid], align, match, mismatch, gap)
25    align ← Hirschberg(X[xmid+1:xlen], Y[ymid+1:ylen], match, mismatch, gap)
26  return align
```
</details>
<details>
<summary>Download files from GitHub</summary>
<details>
<summary>Basic Git settings</summary>

> * Configure the Git editor
> ```bash
> git config --global core.editor notepad
> ```
> * Configure your name and email address
> ```bash
> git config --global user.name "Zuzana Nova"
> git config --global user.email z.nova@vut.cz
> ```
> * Check current settings
> ```bash
> git config --global --list
> ```
>
</details>

* Create a fork on your GitHub account. 
  On the GitHub page of this repository find a <kbd>Fork</kbd> button in the upper right corner.
  
* Clone forked repository from your GitHub page to your computer:
```bash
git clone <fork repository address>
```
* In a local repository, set new remote for a project repository:
```bash
git remote add upstream https://github.com/mpa-prg/exercise_05.git
```

#### Send files to GitHub
Create a new commit and send new changes to your remote repository.
* Add file to a new commit.
```bash
git add <file_name>
```
* Create a new commit, enter commit message, save the file and close it.
```bash
git commit
```
* Send a new commit to your GitHub repository.
```bash
git push origin main
```

</details>
