function [Score, Alignment, Scoringmatrix, Plot] = nw_strottmann(seq1,seq2,AoN,scm,gapopen,gapextend)

if nargin > 0
    seq1 = convertStringsToChars(seq1);
end

if nargin > 1
    seq2 = convertStringsToChars(seq2);
end

AlignmentOptions.scm = scm;
AlignmentOptions.type = AoN;
AlignmentOptions.gapopen = gapopen;
AlignmentOptions.gapextend = gapextend;
AlignmentOptions.lseq1 = numel(seq1);
AlignmentOptions.lseq2 = numel(seq2);

if AlignmentOptions.gapopen == AlignmentOptions.gapextend
    [Score, Alignment, scoring_matrix,x,y] = nw_lineare_gap_kosten(seq1,seq2,AlignmentOptions);
    Scoringmatrix = scoring_matrix;
    S = scoring_matrix;
else
    [Score, Alignment, scoring_matrix,hilfs_matrix_x,hilfs_matrix_y,x,y] = nw_affine_gap_kosten(seq1,seq2,AlignmentOptions);
    S(:,:,1) = scoring_matrix;
    S(:,:,2) = hilfs_matrix_x;
    S(:,:,3) = hilfs_matrix_y;
    Scoringmatrix.scores = S;
end

if nargout > 3
Plot = figure;
abbildungs_matrix = max(S(2:end,2:end,:),[],3);
clim=max(max(max(abs(abbildungs_matrix(~isinf(abbildungs_matrix))))),eps);
imagesc(abbildungs_matrix,[-clim clim]);
colormap(turbo);
set(colorbar,'YLim',[min(abbildungs_matrix(:)) max(abbildungs_matrix(:))])
title('Wertungsraum und optimaler Pfad')
xlabel('Sequenz 1')
ylabel('Sequenz 2')
hold on
plot(x,y,'k.')
end


function [Score, Alignment, scoring_matrix,x,y] = nw_lineare_gap_kosten(seq1,seq2,AlignmentOptions)

%Kosten_Matrix erstellen
kosten_matrix = kosten_matrix_aufstellen(seq1,seq2,AlignmentOptions);

% Scoring- und Richtungsmatrix erstellen
% Richtungen
% 1 = scoring_matrix(i,j) berechnet aus scoring_matrix(i-1,j-1)
% 2 = scoring_matrix(i,j) berechnet aus scoring_matrix(i-1,j)
% 3 = scoring_matrix(i,j) berechnet aus scoring_matrix(i-1,j-1)=(i-1,j)
% 4 = scoring_matrix(i,j) berechnet aus scoring_matrix(i,j-1)
% 5 = scoring_matrix(i,j) berechnet aus scoring_matrix(i-1,j-1)=(i,j-1)

scoring_matrix = zeros(AlignmentOptions.lseq2+1,AlignmentOptions.lseq1+1);
direction_matrix = zeros(AlignmentOptions.lseq2+1,AlignmentOptions.lseq1+1);

for i=2:AlignmentOptions.lseq2+1
    scoring_matrix(i,1) = (i-1)*(-AlignmentOptions.gapopen);
    direction_matrix(i,1) = 2;
end

for j=2:AlignmentOptions.lseq1+1
    scoring_matrix(1,j) = (j-1)*(-AlignmentOptions.gapopen);
    direction_matrix(1,j) = 4;
end

H = [0, 0, 0];

for i=2:AlignmentOptions.lseq2+1
    for j=2:AlignmentOptions.lseq1+1
        H(1) = scoring_matrix(i-1,j-1)+kosten_matrix(i-1,j-1);
        H(2) = scoring_matrix(i-1,j)-AlignmentOptions.gapopen;
        H(3) = scoring_matrix(i,j-1)-AlignmentOptions.gapopen;
        scoring_matrix(i,j) = max(H);
        
        if H(1) == H(2) && H(1) == max(H)
            direction_matrix(i,j) = 3;
        elseif H(1) == H(3) && H(1) == max(H)
            direction_matrix(i,j) = 5;
        elseif H(2)== max(H)
            direction_matrix(i,j) = 2;
        elseif H(3) == max(H)
            direction_matrix(i,j) = 4;
        else
            direction_matrix(i,j) = 1;
        end
    end
end

%Traceback

i = AlignmentOptions.lseq2+1;
j = AlignmentOptions.lseq1+1;
n = 0;
nn = 0;

while direction_matrix(i,j) ~= 0
    if direction_matrix(i,j)==2
        i = i-1;
        n = n+1;
    elseif direction_matrix(i,j)==4
        j = j-1;
        n = n+1;
    else
        i = i-1;
        j = j-1;
        n = n+1;
        nn = nn+1;
    end
end

i = AlignmentOptions.lseq2+1;
j = AlignmentOptions.lseq1+1;
k = 0;
v = 0;
h = 0;
x = zeros(nn);
y = zeros(nn);
SEQ1 = '';
SEQ2 = '';
VERB = '';

while direction_matrix(i,j) ~= 0
    if direction_matrix(i,j)==2
        i = i-1;
        SEQ1(n-h) = '-';
        SEQ2(n-h) = seq2(end-k);
        VERB(n-h) = ' ';
        k = k+1;
        h = h+1;
    elseif direction_matrix(i,j)==4
        j = j-1;
        SEQ1(n-h) = seq1(end-v);
        SEQ2(n-h) = '-';
        VERB(n-h) = ' ';
        v = v+1;
        h = h+1;
    else
        i = i-1;
        j = j-1;
        SEQ1(n-h) = seq1(end-v);
        SEQ2(n-h) = seq2(end-k);
        x(nn) = AlignmentOptions.lseq1-v;
        y(nn) = AlignmentOptions.lseq2-k;
        nn = nn-1;
        if AlignmentOptions.type == 'AA'
            value = AlignmentOptions.scm(bst_zu_zahl_A(SEQ2(n-h)),bst_zu_zahl_A(SEQ1(n-h)));
        else
            value = AlignmentOptions.scm(bst_zu_zahl_N(SEQ2(n-h)),bst_zu_zahl_N(SEQ1(n-h)));
        end
        if SEQ1(n-h) == SEQ2(n-h)
            VERB(n-h) = '|';
        elseif value >= 0
            VERB(n-h) = ':';
        else
            VERB(n-h) = ' ';
        end
        k = k+1;
        v = v+1;
        h = h+1;
    end
end

Alignment = [SEQ1;VERB;SEQ2];
Score = scoring_matrix(AlignmentOptions.lseq2+1,AlignmentOptions.lseq1+1);


function [Score, Alignment, scoring_matrix,hilfs_matrix_x,hilfs_matrix_y,x,y] = nw_affine_gap_kosten(seq1,seq2,AlignmentOptions)

%Kosten_Matrix erstellen
kosten_matrix = kosten_matrix_aufstellen(seq1,seq2,AlignmentOptions);

% Scoring- und Richtungsmatrix erstellen
scoring_matrix = zeros(AlignmentOptions.lseq2+1,AlignmentOptions.lseq1+1);
direction_matrix = zeros(AlignmentOptions.lseq2+1,AlignmentOptions.lseq1+1);
direction_matrix_x = zeros(AlignmentOptions.lseq2+1,AlignmentOptions.lseq1+1);
direction_matrix_y = zeros(AlignmentOptions.lseq2+1,AlignmentOptions.lseq1+1);
hilfs_matrix_x = zeros(AlignmentOptions.lseq2+1,AlignmentOptions.lseq1+1);
hilfs_matrix_y = zeros(AlignmentOptions.lseq2+1,AlignmentOptions.lseq1+1);

scoring_matrix(1,1) = 0;
hilfs_matrix_y(1,1) = -inf;
hilfs_matrix_x(1,1) = -inf;
    
for i=2:AlignmentOptions.lseq2+1
    scoring_matrix(i,1) = -inf;
    hilfs_matrix_y(i,1) = -inf;
    hilfs_matrix_x(i,1) = -AlignmentOptions.gapopen-((i-1)-1)*AlignmentOptions.gapextend;
    direction_matrix_x(i,1) = 2;
end

for j=2:AlignmentOptions.lseq1+1
    scoring_matrix(1,j) = -inf;
    hilfs_matrix_x(1,j) = -inf;
    hilfs_matrix_y(1,j) = -AlignmentOptions.gapopen-((j-1)-1)*AlignmentOptions.gapextend;
    direction_matrix_y(1,j) = 4;
end

Hs = [0, 0, 0];
Hhx = [0, 0];
Hhy = [0, 0];

% Richtungen
% 1 = scoring_matrix/hilfs_matrix_x/hilfs_matrix_y Eintrag berechnet aus 
% vorherigem scoring_matrix Eintrag
% 2 = scoring_matrix/hilfs_matrix_x Eintrag berechnet aus vorherigem
% hilfs_matrix_x Eintrag
% 3 = scoring_matrix/hilfs_matrix_x Eintrag berechnet aus vorherigem 
% scoring_matrix/hilfs_matrix_x Eintrag (M=Ix) 
% 4 = scoring_matrix/hilfs_matrix_y Eintrag berechnet aus vorherigem
% hilfs_matrix_y Eitrag
% 5 = scoring_matrix/hilfs_matrix_y Eintrag berechnet aus vorherigem 
% scoring_matrix/hilfs_matrix_y Eintrag (M=Iy)

for i=2:AlignmentOptions.lseq2+1
    for j=2:AlignmentOptions.lseq1+1
        Hs(1) = scoring_matrix(i-1,j-1)+kosten_matrix(i-1,j-1);
        Hs(2) = hilfs_matrix_x(i-1,j-1)+kosten_matrix(i-1,j-1);
        Hs(3) = hilfs_matrix_y(i-1,j-1)+kosten_matrix(i-1,j-1);
        scoring_matrix(i,j) = max(Hs);
        
        Hhx(1) = scoring_matrix(i-1,j)-AlignmentOptions.gapopen;
        Hhx(2) = hilfs_matrix_x(i-1,j)-AlignmentOptions.gapextend;
        hilfs_matrix_x(i,j) = max(Hhx);
        
        Hhy(1) = scoring_matrix(i,j-1)-AlignmentOptions.gapopen;
        Hhy(2) = hilfs_matrix_y(i,j-1)-AlignmentOptions.gapextend;
        hilfs_matrix_y(i,j) = max(Hhy);
        
        if Hs(1) == Hs(2) && Hs(1) == max(Hs)
            direction_matrix(i,j) = 3;
        elseif Hs(1) == Hs(3) && Hs(1) == max(Hs)
            direction_matrix(i,j) = 5;
        elseif Hs(1) == max(Hs)
             direction_matrix(i,j) = 1;
        elseif(Hs(2)== max(Hs))
             direction_matrix(i,j) = 2;
        else
             direction_matrix(i,j) = 4;
        end
        
        if Hhx(1) == Hhx(2)
             direction_matrix_x(i,j) = 3;
        elseif Hhx(1) == max(Hhx)
             direction_matrix_x(i,j) = 1;
        else
             direction_matrix_x(i,j) = 2;
        end
        
        if Hhy(1) == Hhy(2)
             direction_matrix_y(i,j) = 5;
        elseif Hhy(1) == max(Hhy)
             direction_matrix_y(i,j) = 1;
        else
             direction_matrix_y(i,j) = 4;
        end
        
    end
end

%Traceback

i = AlignmentOptions.lseq2+1;
j = AlignmentOptions.lseq1+1;
n = 0;
nn = 0;
t = 0;
traceback_hilfe = [0,0,0];

traceback_hilfe(1) = scoring_matrix(i,j);
traceback_hilfe(2) = hilfs_matrix_x(i,j);
traceback_hilfe(3) = hilfs_matrix_y(i,j);
    
if traceback_hilfe(1) == max(traceback_hilfe)
    tracebackmatrix = direction_matrix;
    t = 1;
elseif traceback_hilfe(2) == max(traceback_hilfe)
    tracebackmatrix = direction_matrix_x;
    t = 2;
else
    tracebackmatrix = direction_matrix_y;
    t = 3;
end

while i ~= 1 || j ~= 1
    if tracebackmatrix(i,j) == 1 && t == 1
        j = j-1;
        i = i-1; 
        n = n+1;
        nn = nn+1; 
    elseif tracebackmatrix(i,j) == 3 && t == 1
        j = j-1;
        i = i-1; 
        n = n+1;
        nn = nn+1;
        tracebackmatrix = direction_matrix_x;
        t = 2;
    elseif tracebackmatrix(i,j) == 5 && t == 1
        j = j-1;
        i = i-1; 
        n = n+1;
        nn = nn+1;
        tracebackmatrix = direction_matrix_y;
        t = 3;
    elseif tracebackmatrix(i,j) == 2 && t == 1
        j = j-1;
        i = i-1; 
        n = n+1;
        nn = nn+1;
        tracebackmatrix = direction_matrix_x;
        t = 2;
    elseif tracebackmatrix(i,j) == 4 && t == 1
        j = j-1;
        i = i-1; 
        n = n+1;
        nn = nn+1;
        tracebackmatrix = direction_matrix_y;
        t = 3;
    elseif tracebackmatrix(i,j) == 1 && t == 2
        i = i-1;
        n = n+1;
        tracebackmatrix = direction_matrix;
        t = 1;
    elseif tracebackmatrix(i,j) == 2 && t == 2
        i = i-1;
        n = n+1;
    elseif tracebackmatrix(i,j) == 3 && t == 2
        i = i-1;
        n = n+1;
    elseif tracebackmatrix(i,j) == 1 && t == 3
        j = j-1;
        n = n+1;
        tracebackmatrix = direction_matrix;
        t = 1;
    elseif tracebackmatrix(i,j) == 5 && t == 3
        j = j-1;
        n = n+1;
    elseif tracebackmatrix(i,j) == 4 && t == 3
        j = j-1;
        n = n+1;
    end
end

i = AlignmentOptions.lseq2+1;
j = AlignmentOptions.lseq1+1;
k = 0;
v = 0;
h = 0;
t = 0;
x = zeros(nn);
y = zeros(nn);
SEQ1 = '';
SEQ2 = '';
VERB = '';
    
if traceback_hilfe(1) == max(traceback_hilfe)
    tracebackmatrix = direction_matrix;
    t = 1;
elseif traceback_hilfe(2) == max(traceback_hilfe)
    tracebackmatrix = direction_matrix_x;
    t = 2;
else
    tracebackmatrix = direction_matrix_y;
    t = 3;
end

while i ~= 1 || j ~= 1
    if tracebackmatrix(i,j) == 1 && t == 1
        j = j-1;
        i = i-1; 
        SEQ1(n-h) = seq1(end-v);
        SEQ2(n-h) = seq2(end-k);
        x(nn) = AlignmentOptions.lseq1-v;
        y(nn) = AlignmentOptions.lseq2-k;
        nn = nn-1;
        if AlignmentOptions.type == 'AA'
            value = AlignmentOptions.scm(bst_zu_zahl_A(SEQ2(n-h)),bst_zu_zahl_A(SEQ1(n-h)));
        else
            value = AlignmentOptions.scm(bst_zu_zahl_N(SEQ2(n-h)),bst_zu_zahl_N(SEQ1(n-h)));
        end
        if SEQ1(n-h) == SEQ2(n-h)
            VERB(n-h) = '|';
        elseif value >= 0
            VERB(n-h) = ':';
        else
            VERB(n-h) = ' ';
        end
        k = k+1;
        v = v+1;
        h = h+1;
    elseif tracebackmatrix(i,j) == 3 && t == 1
        j = j-1;
        i = i-1; 
        SEQ1(n-h) = seq1(end-v);
        SEQ2(n-h) = seq2(end-k);
        x(nn) = AlignmentOptions.lseq1-v;
        y(nn) = AlignmentOptions.lseq2-k;
        nn = nn-1;
        if AlignmentOptions.type == 'AA'
            value = AlignmentOptions.scm(bst_zu_zahl_A(SEQ2(n-h)),bst_zu_zahl_A(SEQ1(n-h)));
        else
            value = AlignmentOptions.scm(bst_zu_zahl_N(SEQ2(n-h)),bst_zu_zahl_N(SEQ1(n-h)));
        end
        if SEQ1(n-h) == SEQ2(n-h)
            VERB(n-h) = '|';
        elseif value >= 0
            VERB(n-h) = ':';
        else
            VERB(n-h) = ' ';
        end
        k = k+1;
        v = v+1;
        h = h+1;
        tracebackmatrix = direction_matrix_x;
        t = 2;
    elseif tracebackmatrix(i,j) == 5 && t == 1
        j = j-1;
        i = i-1; 
        SEQ1(n-h) = seq1(end-v);
        SEQ2(n-h) = seq2(end-k);
        x(nn) = AlignmentOptions.lseq1-v;
        y(nn) = AlignmentOptions.lseq2-k;
        nn = nn-1;
        if AlignmentOptions.type == 'AA'
            value = AlignmentOptions.scm(bst_zu_zahl_A(SEQ2(n-h)),bst_zu_zahl_A(SEQ1(n-h)));
        else
            value = AlignmentOptions.scm(bst_zu_zahl_N(SEQ2(n-h)),bst_zu_zahl_N(SEQ1(n-h)));
        end
        if SEQ1(n-h) == SEQ2(n-h)
            VERB(n-h) = '|';
        elseif value >= 0
            VERB(n-h) = ':';
        else
            VERB(n-h) = ' ';
        end
        k = k+1;
        v = v+1;
        h = h+1;
        tracebackmatrix = direction_matrix_y;
        t = 3;
    elseif tracebackmatrix(i,j) == 2 && t == 1
        j = j-1;
        i = i-1; 
        SEQ1(n-h) = seq1(end-v);
        SEQ2(n-h) = seq2(end-k);
        x(nn) = AlignmentOptions.lseq1-v;
        y(nn) = AlignmentOptions.lseq2-k;
        nn = nn-1;
        if AlignmentOptions.type == 'AA'
            value = AlignmentOptions.scm(bst_zu_zahl_A(SEQ2(n-h)),bst_zu_zahl_A(SEQ1(n-h)));
        else
            value = AlignmentOptions.scm(bst_zu_zahl_N(SEQ2(n-h)),bst_zu_zahl_N(SEQ1(n-h)));
        end
        if SEQ1(n-h) == SEQ2(n-h)
            VERB(n-h) = '|';
        elseif value >= 0
            VERB(n-h) = ':';
        else
            VERB(n-h) = ' ';
        end
        k = k+1;
        v = v+1;
        h = h+1;
        tracebackmatrix = direction_matrix_x;
        t = 2;
    elseif tracebackmatrix(i,j) == 4 && t == 1
        j = j-1;
        i = i-1; 
        SEQ1(n-h) = seq1(end-v);
        SEQ2(n-h) = seq2(end-k);
        x(nn) = AlignmentOptions.lseq1-v;
        y(nn) = AlignmentOptions.lseq2-k;
        nn = nn-1;
        if AlignmentOptions.type == 'AA'
            value = AlignmentOptions.scm(bst_zu_zahl_A(SEQ2(n-h)),bst_zu_zahl_A(SEQ1(n-h)));
        else
            value = AlignmentOptions.scm(bst_zu_zahl_N(SEQ2(n-h)),bst_zu_zahl_N(SEQ1(n-h)));
        end
        if SEQ1(n-h) == SEQ2(n-h)
            VERB(n-h) = '|';
        elseif value >= 0
            VERB(n-h) = ':';
        else
            VERB(n-h) = ' ';
        end
        k = k+1;
        v = v+1;
        h = h+1;
        tracebackmatrix = direction_matrix_y;
        t = 3;
    elseif tracebackmatrix(i,j) == 1 && t == 2
        i = i-1;
        SEQ1(n-h) = '-';
        SEQ2(n-h) = seq2(end-k);
        VERB(n-h) = ' ';
        k = k+1;
        h = h+1;
        tracebackmatrix = direction_matrix;
        t = 1;
    elseif tracebackmatrix(i,j) == 2 && t == 2
        i = i-1;
        SEQ1(n-h) = '-';
        SEQ2(n-h) = seq2(end-k);
        VERB(n-h) = ' ';
        k = k+1;
        h = h+1;
    elseif tracebackmatrix(i,j) == 3 && t == 2
        i = i-1;
        SEQ1(n-h) = '-';
        SEQ2(n-h) = seq2(end-k);
        VERB(n-h) = ' ';
        k = k+1;
        h = h+1;
    elseif tracebackmatrix(i,j) == 1 && t == 3
        j = j-1;
        SEQ1(n-h) = seq1(end-v);
        SEQ2(n-h) = '-';
        VERB(n-h) = ' ';
        v = v+1;
        h = h+1;
        tracebackmatrix = direction_matrix;
        t = 1;
    elseif tracebackmatrix(i,j) == 5 && t == 3
        j = j-1;
        SEQ1(n-h) = seq1(end-v);
        SEQ2(n-h) = '-';
        VERB(n-h) = ' ';
        v = v+1;
        h = h+1;
    elseif tracebackmatrix(i,j) == 4 && t == 3
        j = j-1;
        SEQ1(n-h) = seq1(end-v);
        SEQ2(n-h) = '-';
        VERB(n-h) = ' ';
        v = v+1;
        h = h+1;
    end
end

Alignment = [SEQ1;VERB;SEQ2];
score_hilfe = [0,0,0];
score_hilfe(1) = scoring_matrix(AlignmentOptions.lseq2+1,AlignmentOptions.lseq1+1);
score_hilfe(2) = hilfs_matrix_x(AlignmentOptions.lseq2+1,AlignmentOptions.lseq1+1);
score_hilfe(3) = hilfs_matrix_y(AlignmentOptions.lseq2+1,AlignmentOptions.lseq1+1);
Score = max(score_hilfe);


function A = bst_zu_zahl_A(sequenz)
A = zeros(1,numel(sequenz));

for i = 1:numel(sequenz)
    if sequenz(i) == 'A'
        A(i) = 1;
    elseif sequenz(i) == 'R'
        A(i) = 2;
    elseif sequenz(i) == 'N'
        A(i) = 3;
    elseif sequenz(i) == 'D'
        A(i) = 4;
    elseif sequenz(i) == 'C'
        A(i) = 5;
    elseif sequenz(i) == 'Q'
        A(i) = 6;
    elseif sequenz(i) == 'E'
        A(i) = 7;
    elseif sequenz(i) == 'G'
        A(i) = 8;
    elseif sequenz(i) == 'H'
        A(i) = 9;
    elseif sequenz(i) == 'I'
        A(i) = 10;
    elseif sequenz(i) == 'L'
        A(i) = 11;
    elseif sequenz(i) == 'K'
        A(i) = 12;
    elseif sequenz(i) == 'M'
        A(i) = 13;
    elseif sequenz(i) == 'F'
        A(i) = 14;
    elseif sequenz(i) == 'P'
        A(i) = 15;
    elseif sequenz(i) == 'S'
        A(i) = 16;
    elseif sequenz(i) == 'T'
        A(i) = 17;
    elseif sequenz(i) == 'W'
        A(i) = 18;
    elseif sequenz(i) == 'Y'
        A(i) = 19;
    elseif sequenz(i) == 'V'
        A(i) = 20;
    elseif sequenz(i) == 'B'
        A(i) = 21;
    elseif sequenz(i) == 'Z'
        A(i) = 22;
    elseif sequenz(i) == 'X'
        A(i) = 23;
    elseif sequenz(i) == '*'
        A(i) = 24;
    else 
        A(i) = 0;
    end
end


function A = bst_zu_zahl_N(sequenz)
A = zeros(1,numel(sequenz));

for i = 1:numel(sequenz)
    if sequenz(i) == 'A'
        A(i) = 1;
    elseif sequenz(i) == 'C'
        A(i) = 2;
    elseif sequenz(i) == 'G'
        A(i) = 3;
    elseif sequenz(i) == 'T'
        A(i) = 4;
    elseif sequenz(i) == 'R'
        A(i) = 5;
    elseif sequenz(i) == 'Y'
        A(i) = 6;
    elseif sequenz(i) == 'K'
        A(i) = 7;
    elseif sequenz(i) == 'M'
        A(i) = 8;
    elseif sequenz(i) == 'S'
        A(i) = 9;
    elseif sequenz(i) == 'W'
        A(i) = 10;
    elseif sequenz(i) == 'B'
        A(i) = 11;
    elseif sequenz(i) == 'D'
        A(i) = 12;
    elseif sequenz(i) == 'H'
        A(i) = 13;
    elseif sequenz(i) == 'V'
        A(i) = 14;
    elseif sequenz(i) == 'N'
        A(i) = 15;
    else 
        A(i) = 0;
    end
end


function k_matrix = kosten_matrix_aufstellen(s1,s2,AlignmentOptions)
if AlignmentOptions.type == 'AA'
    B = bst_zu_zahl_A(s1);
    C = bst_zu_zahl_A(s2);
else
    B = bst_zu_zahl_N(s1);
    C = bst_zu_zahl_N(s2);
end
    k_matrix = zeros(numel(s2),numel(s1));

for i = 1:numel(s2)
    for j = 1: numel(s1)
        k_matrix(i,j) = AlignmentOptions.scm(C(i),B(j));
    end
end
