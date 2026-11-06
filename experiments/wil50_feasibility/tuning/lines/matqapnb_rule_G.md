# MatQapNB's branching rule G (the Gilmore-Lawler rule), as coded

Source: `NewtBracket-v2.0/applications/MatQapNB` (Fujii, Ito, Kim, Kojima,
Shinano, Toh, "Solving Challenging Large Scale QAPs", ZIB Report 21-02),
Nathan's copy in `~/Downloads`. Rule G is `paramNB.branchSW = 4`
(`BandBexamples/QapExample31.m`: "The branching rule G was not mentioned
there but utilizes the Gilmore-Lawler bound of Q((F,f),(L,l)) for phi(f,l)
in the branching rules described in Section 5"). The default is rule P
(`branchSW = 2`, `defaultParamsBrQAP.m:38`). Code below is verbatim from
`branching/branch.m` (case 4 of `genBranchCandidate`), `branching/GLbound.m`
and `subPrograms/genAssignMat.m`, with the dead comment lines removed.

## What it does, in words

At a node with assignment `facToLoc` (facility -> location, 0 = free), with
free facilities F and free locations L (|F| = |L| = q):

1. For every free pair (f, l), form the child node that adds f -> l and
   compute phi(f, l) = the Gilmore-Lawler bound of that child. The GL bound
   is `GLbound`: reduce to the q-1 free facilities and locations, carry the
   fixed-fixed part as a constant, the fixed-free part as a linear term
   (2 * row sums of Q over the fixed pairs), form the GL cost matrix
   `CMat = bHat * aHat' + kron(diag B', diag A) + linear` (rows of A sorted
   ascending, rows of B descending, each without its diagonal: the classical
   minimal scalar product), and solve the LAP on it with `munkres`. That is q^2
   LAPs of size q-1 per branching decision.
2. `objMat(i, j) = phi(F(i), L(j))`. Take row means (per facility) and column
   means (per location). The facility with the largest mean GL bound over its
   locations, or the location with the largest mean over its facilities,
   whichever mean is larger, is the branching line. (Ties go to the facility.)
3. Children: only the pairs on that line whose `ceil(phi) < bestObjVal`.
   The others are pruned at creation (`objActiveLocIdx01` / `objActiveFacIdx01`).
   Rule D (case 3) uses the same mean-then-prune shape with the dual DNN
   bound as phi; rule M (case 1) and rule P (case 2) take the mean/primal
   projection and create all q children of the line without pruning.
4. Each child is then given `LBv = max(parent LBv, GLbound(child))`
   (`selectBranch`), so the GL bound is also every child's initial bound,
   under every rule.

So rule G is "branch on the facility (or location) whose children have the
largest average Gilmore-Lawler bound, and drop the children GL already
fathoms". It is the mean of the child bounds, not the min or a product rule.

## Code

```matlab
% branching/branch.m, genBranchCandidate, case 4 ('GLB')
facToLoc0 = oneNode.facToLoc;
qapDim0 = length(facToLoc0);
[I,~,V] = find(facToLoc0);
availableFacIdx = sort(setdiff((1:qapDim0)',I),'ascend');
availableLocIdx = sort(setdiff((1:qapDim0)',V),'ascend');

objMatSize = length(availableLocIdx);
cellSize = objMatSize*objMatSize;
objVect = zeros(1,cellSize);
facToLocCell = cell(1,cellSize);
k = 0;
for j=1:objMatSize
    q = availableLocIdx(j);
    for i=1:objMatSize
        k = k+1;
        facToLocCell{k} = facToLoc0;
        p = availableFacIdx(i);
        facToLocCell{k}(p) = q;          % child: facility p -> location q
    end
end
for k=1:cellSize
    [objVect(k)] = GLbound(QMat0,A0,B0,facToLocCell{k});
end
objMat = reshape(objVect,objMatSize,objMatSize);   % rows: facilities, cols: locations

objActiveFacColVec01 = zeros(objMatSize,1);
objActiveLocIdx01 = cell(1,objMatSize);
for i=1:objMatSize
    objActiveLocIdx01{i} = find(ceil(objMat(i,:)) < paramNB.bestObjVal);
    objActiveFacColVec01(i) = length(objActiveLocIdx01{i});
end
objActiveLocRowVec01 = zeros(1,objMatSize);
objActiveFacIdx01 = cell(1,objMatSize);
for j=1:objMatSize
    objActiveFacIdx01{j} = find(ceil(objMat(:,j)) < paramNB.bestObjVal);
    objActiveLocRowVec01(j) = length(objActiveFacIdx01{j});
end

objMeanFacColVec = mean(objMat,2);   % mean GL bound of each facility's children
objMeanLocRowVec = mean(objMat,1);   % mean GL bound of each location's children
[maxFacVal,maxFacLidx] = max(objMeanFacColVec(:));
[maxLocVal,maxLocLidx] = max(objMeanLocRowVec(:));
if maxFacVal >= maxLocVal
    chosenFacLidxSet = maxFacLidx;                           % branch on this facility
    chosenLocLidxSet = objActiveLocIdx01{maxFacLidx}(:)';    % its unfathomed locations
else
    chosenFacLidxSet = objActiveFacIdx01{maxLocLidx}(:)';    % its unfathomed facilities
    chosenLocLidxSet = maxLocLidx;                           % branch on this location
end

if isempty(chosenLocLidxSet), oneNode.locGidxSet = availableLocIdx;
else, oneNode.locGidxSet = availableLocIdx(chosenLocLidxSet); end
if isempty(chosenFacLidxSet), oneNode.facGidxSet = availableFacIdx;
else, oneNode.facGidxSet = availableFacIdx(chosenFacLidxSet); end
% selectBranch then creates one child per (f in facGidxSet, l in locGidxSet)
% with LBv = max(parent.LBv, GLbound(child)).
```

```matlab
% branching/GLbound.m
function [LBv] = GLbound(QMat0,A0,B0,facToLoc)
    qapDim0 = length(facToLoc(:));
    [assignMat] = genAssignMat(facToLoc);
    [I,~,V] = find(facToLoc);
    availableFacIdx = sort(setdiff((1:qapDim0)',I),'ascend');
    availableLocIdx = sort(setdiff((1:qapDim0)',V),'ascend');
    A = A0(availableLocIdx,availableLocIdx);
    B = B0(availableFacIdx,availableFacIdx);
    qapDim = length(availableFacIdx);
    [linearVect,objConst] = subQAP(QMat0,qapDim,assignMat);
    [CMat] = linearlizeQAP(qapDim,A,B,linearVect);
    [~,cost] = munkres(CMat);
    LBv = cost + objConst;
end

function [linearVect,objConst,QMat] = subQAP(QMat0,qapDim,assignMat)
    qapDim2 = qapDim*qapDim;
    I1 = find(assignMat(:) == 1);     % fixed pairs (as indices into the n^2 vector)
    Ia = find(assignMat(:) == 0);     % free pairs
    if ~isempty(I1)
        objConst = sum(sum(QMat0(I1,I1),1));
        linearVect = 2*sum(QMat0(Ia,I1),2);
    else
        objConst = 0;
        linearVect = zeros(qapDim2,1);
    end
    if nargout == 3, QMat = QMat0(Ia,Ia); end
end

function [CMat] = linearlizeQAP(qapDim,A,B,linearVect)
    aHatMat = zeros(qapDim,qapDim-1);
    bHatMat = zeros(qapDim,qapDim-1);
    for i=1:qapDim
        ahat = A(i,[1:(i-1),(i+1):qapDim]);
        aHatMat(i,:) = sort(ahat,'ascend');
    end
    for j=1:qapDim
        bhat = B(j,[1:(j-1),(j+1):qapDim]);
        bHatMat(j,:) = sort(bhat,'descend');
    end
    diagA = diag(A);
    diagB = diag(B);
    CMat = bHatMat*aHatMat' + kron(diagB',diagA) + reshape(linearVect,qapDim,qapDim);
end
```

```matlab
% subPrograms/genAssignMat.m: 1 = fixed pair, -1 = excluded by a fixed row or
% column, 0 = free
function [assignMat] = genAssignMat(facToLoc)
    qapDim = length(facToLoc);
    assignMat = zeros(qapDim,qapDim);
    [I,~,K] = find(facToLoc);
    assignMat(I,:) = -1;
    assignMat(:,K) = -1;
    for p=1:length(I)
        assignMat(I(p),K(p)) = 1;
    end
end
```

`QMat0` is the n^2 x n^2 objective matrix of the QAP in the lifted form the
package uses (`subPrograms/QAPtoDNNformatBr.m`); `A0`, `B0` are the distance
and flow matrices; `munkres` is the bundled Hungarian solver
(`solver/munkres/`). Note `A` is indexed by locations and `B` by facilities in
`GLbound`, and the GL matrix is built as `bHat * aHat'`, so `CMat(i, j)` is the
cost of free facility `availableFacIdx(i)` at free location
`availableLocIdx(j)`; the rows of `objMat` run over facilities, the columns
over locations, consistent with that.

## Relation to our rules

Our rule 6 (eccentric), rule 4/7 and the lookahead rule 9 pick a line from
the GL or QPB child bounds too, but by the strongest-children criterion
(the min or the number of fathomed children) rather than the mean; rule 8
reads the DNN marginals. Rule G's mean-over-children criterion with GL
pruning is cheap (q^2 LAPs of size q-1, about 2,000 LAPs of size 45 at
wil50 depth 4) and is one of the candidates the Nov 6-7 tail line study
(`lines.py --rules M,D`) can take next if M or D do anything.
