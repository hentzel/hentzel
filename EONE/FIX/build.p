
program build(input,output);

const PRIME = 251; (* NOT USED *)

      c1 = 2;
      c2 = 3;
      c3 = 7;
      c4 = 7;

type poly = array[0..1000,1..2] of longint;       (* coefficient basis     *)
                                                   (*  p[0,1] is the length *)
var
    i,j,k,n:longint;    
    H:array[0..140000000,1..2] of longint;
             (*  139799736   *)
    index:array[0.. 503721 ,1..4] of longint;
    name,filein,fileout:text;
    A,B,C,D,E,F,G,Z,ans:poly;     
    AA:array[0..100000] of poly;

    mat:array[1..14502] of longint;
    col:longint;
    i1,i2,i3,i4,i5,i6,i7,i8:longint;
    x1,x2,x3,x4,x5,x6,x7,x8:poly;
    w:array[1..8] of longint;
    start,stop:array[1..16] of longint;
    comstart,comstop:array[1..16] of longint;
    
    T:array[1..55653] of string[60];
    xreal:double;
    degone,degtwo:longint;

    combasis:array[1..55653,1..4] of longint;
    basis:array[1..70000,1..4] of longint;
    log:array[1..70000] of longint;
    comlog:array[1..55653] of longint;
    combasislength:longint;
    eonebasislength:longint;
    coefficient:array[1..55653] of longint;
      etermfile:array[1..55653] of longint;
       comparts:array[1..55653,1..5] of longint; 
           comT:array[1..55653] of string[60];

function degree(n:longint):longint;

var i,deg:longint;

begin
    deg := 0;
    for i:=1 to 8 do
      if (start[i] <= n)and(n<=stop[i]) then
          deg := i;
    if deg = 0  then
       begin
         writeln(' PROBLEM, Degree is zero ');
         writeln(' n is ',n,'    PAUSED TILL THE COWS COME HOME ');
         readln(deg);
       end;
    degree := deg;
end;

procedure displaypoly(A:poly);

var i:longint;

begin
   writeln(' The length of the poly is ',A[0,1]:10);
   for i:=1 to A[0,1] do
      begin
          writeln(A[i,1]:3,' ',A[i,2]:3);
      end; 
end;


procedure seepoly(A:poly);

var i:longint;

begin
   writeln(' The length of the poly is ',A[0,1]:10);
   for i:=1 to A[0,1] do
      begin
          writeln(A[i,1]:3,' ',T[A[i,2]]);
      end; 
end;



function simplify(C:poly):poly;

var i,j:longint;
    D:poly;

begin
    for i:=1 to C[0,1]-1 do
      for j:=i+1 to C[0,1] do
        if C[i,2] = C[j,2] then
           begin
             C[i,1] := (C[i,1]+C[j,1]) mod PRIME;
             C[j,1] := 0;
           end;

    D[0,1] := 0;
    for i:=1 to C[0,1] do
      if C[i,1] <> 0 then
        begin
          D[0,1] := D[0,1]+1;
          D[ D[0,1] ] := C[i];
        end; 
    simplify := D;
end;


(* concats the files and simplies them *)

function xApyB(x:longint;A:poly;y:longint;B:poly):poly;

var i:longint;
    C:poly;

begin
   C[0,1] := 0;
   for i:=1 to A[0,1] do
      begin
        C[0,1] := C[0,1]+1;
        C[ C[0,1], 1 ] := x*A[i,1];
        C[ C[0,1], 2 ] := A[i,2];
      end;
   for i:=1 to B[0,1] do
      begin
        C[0,1] := C[0,1]+1;
        C[ C[0,1], 1 ] := y*B[i,1];
        C[ C[0,1], 2 ] := B[i,2];
      end;
    xApyB := C;
  xApyB := simplify(C);
end;



(* Index    x * y = start to stop                       *)
(* The start to stop terms are in H                     *)
(* This assumes the coefficients of x and y are both 1  *)


function find(x,y:longint):poly;

var n,i:longint;
    A:poly;

begin
    A[0,1] := 0;
    for n:=1 to index[0,1] do
      begin
        if (index[n,1] = x) and (index[n,2] = y) then
           begin
             for i := index[n,3] to index[n,4] do
                begin
                  A[0,1] := A[0,1]+1;
                  A[ A[0,1],1 ] := H[i,1];
                  A[ A[0,1],2 ] := H[i,2];
                end;
           end;
      end;
    find := A;
end;


(* This multiplies polynomials and includes the coefficients *)   
  
function xx(A,B:poly):poly;

var C,ans:poly;
    i,j,k:longint;

begin
   ans[0,1] := 0;
   for i:=1 to A[0,1] do
      for j := 1 to B[0,1] do
         begin
           C := find( A[i,2],B[j,2]);
           for k:=1 to C[0,1] do
              begin
                ans[0,1] := ans[0,1]+1;
                ans[ ans[0,1] ,1] := (C[k,1]*A[i,1]*B[j,1]) mod PRIME; 
                ans[ ans[0,1] ,2] := C[k,2]; 
              end;
            ans := simplify(ans);
          end;
   xx :=  ans;
end;     

function bd(A,B:poly):poly;

begin
    bd := xApyB(1,xx(A,B),-1,xx(B,A));
end;


function pq(A,B,C:poly):poly;

var abxc,axbc:poly;

begin
abxc := xx(xx(A,B),C);
axbc := xx(A,xx(B,C));
abxc := simplify(abxc);
axbc := simplify(axbc);

    pq := xApyB(1,abxc,-1,axbc);
end;



function oo(A,B:poly):poly;

begin
    oo := xApyB(+1,xx(A,B),+1,xx(B,A));
end;


(*   (ab+ba)c + c(ab+ba) - a(bc+cb) - (bc+cb)a     *)

function poq(A,B,C:poly):poly;

begin
       poq :=  xApyB(1, oo(oo(A,B),C),-1,oo(A,oo(B,C)));
end;


function eone(A,B,C,D:poly):poly;

begin
       eone := xApyB(1,xx(pq(A,B,C),D),1,xApyB(1,xx(pq(C,B,D),A),1,xx(pq(D,B,A),C)));
end;

function jone(A,B,C,D:poly):poly;

begin
       jone := xApyB(1,pq(xx(A,B),C,D),1,xApyB(1,pq(xx(B,D),C,A),1,pq(xx(D,A),C,B)));
end;

function xxxx(A,B,C,D:poly):poly;
var ans:poly;

begin
      ans :=               xx(xx(xx(A,B),C),D);
      ans := xApyB(1,ans,1,xx(xx(xx(A,B),D),C));
      ans := xApyB(1,ans,1,xx(xx(xx(A,C),B),D));
      ans := xApyB(1,ans,1,xx(xx(xx(A,C),D),B));
      ans := xApyB(1,ans,1,xx(xx(xx(A,D),B),C));
      ans := xApyB(1,ans,1,xx(xx(xx(A,D),C),B));

      ans := xApyB(1,ans,1,xx(xx(xx(B,C),A),D));
      ans := xApyB(1,ans,1,xx(xx(xx(B,C),D),A));
      ans := xApyB(1,ans,1,xx(xx(xx(B,D),A),C));
      ans := xApyB(1,ans,1,xx(xx(xx(B,D),C),A));

      ans := xApyB(1,ans,1,xx(xx(xx(C,D),A),B));
      ans := xApyB(1,ans,1,xx(xx(xx(C,D),B),A));
      xxxx := ans;
end;
(*---------------------------------------------------------------------*)
procedure xxxxabcd(d1,d2,d3,d4:longint);

var A,B,C,D,ans:poly;
    i1,i2,i3,i4:longint;

begin
     A[0,1] := 1; B[0,1] := 1; C[0,1] := 1; D[0,1] := 1;
     A[1,1] := 1; B[1,1] := 1; C[1,1] := 1; D[1,1] := 1;
     for i1 := start[d1] to stop[d1] do
     for i2 := start[d2] to stop[d2] do
     for i3 := start[d3] to stop[d3] do
     for i4 := start[d4] to stop[d4] do
       begin
          A[1,2] := i1;
          B[1,2] := i2;
          C[1,2] := i3;
          D[1,2] := i4;
          if (d1<>d2)or((d1=d2)and(i1<=i2)) then
          if (d2<>d3)or((d2=d3)and(i2<=i3)) then
          if (d3<>d4)or((d3=d4)and(i3<=i4)) then
            begin
              ans := xxxx(A,B,C,D);
              if ans[0,1] <> 0 then
                begin
                  AA[0][0,1] := AA[0][0,1]+1;
                  AA[ AA[0][0,1] ] := ans;
writeln(name,'"xxxx(',T[i1],',',T[i2],',',T[i3],',',T[i4],')",');
                end;
            end;
       end;         
end;


procedure xxxxabcdx(d1,d2,d3,d4,d5:longint);

var A,B,C,D,E,ans:poly;
    i1,i2,i3,i4,i5:longint;

begin
     A[0,1] := 1; B[0,1] := 1; C[0,1] := 1; D[0,1] := 1; E[0,1] := 1;
     A[1,1] := 1; B[1,1] := 1; C[1,1] := 1; D[1,1] := 1; E[1,1] := 1;
     for i1 := start[d1] to stop[d1] do
     for i2 := start[d2] to stop[d2] do
     for i3 := start[d3] to stop[d3] do
     for i4 := start[d4] to stop[d4] do
     for i5 := start[d5] to stop[d5] do
       begin
          A[1,2] := i1;
          B[1,2] := i2;
          C[1,2] := i3;
          D[1,2] := i4;
          E[1,2] := i5;
          if (d1<>d2)or((d1=d2)and(i1<=i2)) then
          if (d2<>d3)or((d2=d3)and(i2<=i3)) then
          if (d3<>d4)or((d3=d4)and(i3<=i4)) then
            begin
              ans := xx(xxxx(A,B,C,D),E);
              if ans[0,1] <> 0 then
                begin
                  AA[0][0,1] := AA[0][0,1]+1;
                  AA[ AA[0][0,1] ] := ans;
writeln(name,'"(xxxx(',T[i1],',',T[i2],',',T[i3],',',T[i4],')',T[i5],')",');
                end;
            end;
       end;         
end;

procedure xxxxabcdxy(d1,d2,d3,d4,d5,d6:longint);

var A,B,C,D,E,F,ans:poly;
    i1,i2,i3,i4,i5,i6:longint;

begin
     A[0,1] := 1; B[0,1] := 1; C[0,1] := 1; D[0,1] := 1; E[0,1] := 1; F[0,1] := 1;
     A[1,1] := 1; B[1,1] := 1; C[1,1] := 1; D[1,1] := 1; E[1,1] := 1; F[1,1] := 1;
     for i1 := start[d1] to stop[d1] do
     for i2 := start[d2] to stop[d2] do
     for i3 := start[d3] to stop[d3] do
     for i4 := start[d4] to stop[d4] do
     for i5 := start[d5] to stop[d5] do
     for i6 := start[d6] to stop[d6] do
       begin
          A[1,2] := i1;
          B[1,2] := i2;
          C[1,2] := i3;
          D[1,2] := i4;
          E[1,2] := i5;
          F[1,2] := i6;
          if (d1<>d2)or((d1=d2)and(i1<=i2)) then
          if (d2<>d3)or((d2=d3)and(i2<=i3)) then
          if (d3<>d4)or((d3=d4)and(i3<=i4)) then
            begin
              ans := xx(xx(xxxx(A,B,C,D),E),F);
              if ans[0,1] <> 0 then
                begin
                  AA[0][0,1] := AA[0][0,1]+1;
                  AA[ AA[0][0,1] ] := ans;
writeln(name,'"((xxxx(',T[i1],',',T[i2],',',T[i3],',',T[i4],')',T[i5],')',T[i6],')",');
                end;
            end;
       end;         
end;


procedure xxxxabcdxyz(d1,d2,d3,d4,d5,d6,d7:longint);

var A,B,C,D,E,F,G,ans:poly;
    i1,i2,i3,i4,i5,i6,i7:longint;

begin
     A[0,1]:=1;B[0,1]:=1;C[0,1]:=1;D[0,1]:=1;E[0,1]:=1;F[0,1]:=1;G[0,1]:=1;
     A[1,1]:=1;B[1,1]:=1;C[1,1]:=1;D[1,1]:=1;E[1,1]:=1;F[1,1]:=1;G[1,1]:=1;
     for i1 := start[d1] to stop[d1] do
     for i2 := start[d2] to stop[d2] do
     for i3 := start[d3] to stop[d3] do
     for i4 := start[d4] to stop[d4] do
     for i5 := start[d5] to stop[d5] do
     for i6 := start[d6] to stop[d6] do
     for i7 := start[d7] to stop[d7] do
       begin
          A[1,2] := i1;
          B[1,2] := i2;
          C[1,2] := i3;
          D[1,2] := i4;
          E[1,2] := i5;
          F[1,2] := i6;
          G[1,2] := i7;
          if (d1<>d2)or((d1=d2)and(i1<=i2)) then
          if (d2<>d3)or((d2=d3)and(i2<=i3)) then
          if (d3<>d4)or((d3=d4)and(i3<=i4)) then
            begin
              ans := xx(xx(xx(xxxx(A,B,C,D),E),F),G);
              if ans[0,1] <> 0 then
                begin
                  AA[0][0,1] := AA[0][0,1]+1;
                  AA[ AA[0][0,1] ] := ans;
writeln(name,'"(((xxxx(',T[i1],',',T[i2],',',T[i3],',',T[i4],')',T[i5],')',T[i6],')',T[i7],')",');
                end;
            end;
       end;         
end;

procedure xxxxabcdxyzw(d1,d2,d3,d4,d5,d6,d7,d8:longint);

var A,B,C,D,E,F,G,Z,ans:poly;
    i1,i2,i3,i4,i5,i6,i7,i8:longint;

begin
     A[0,1]:=1;B[0,1]:=1;C[0,1]:=1;D[0,1]:=1;E[0,1]:=1;F[0,1]:=1;G[0,1]:=1;Z[0,1]:=1;
     A[1,1]:=1;B[1,1]:=1;C[1,1]:=1;D[1,1]:=1;E[1,1]:=1;F[1,1]:=1;G[1,1]:=1;Z[1,1]:=1;
     for i1 := start[d1] to stop[d1] do
     for i2 := start[d2] to stop[d2] do
     for i3 := start[d3] to stop[d3] do
     for i4 := start[d4] to stop[d4] do
     for i5 := start[d5] to stop[d5] do
     for i6 := start[d6] to stop[d6] do
     for i7 := start[d7] to stop[d7] do
     for i8 := start[d8] to stop[d8] do
       begin
          A[1,2] := i1;
          B[1,2] := i2;
          C[1,2] := i3;
          D[1,2] := i4;
          E[1,2] := i5;
          F[1,2] := i6;
          G[1,2] := i7;
          Z[1,2] := i8;
          if (d1<>d2)or((d1=d2)and(i1<=i2)) then
          if (d2<>d3)or((d2=d3)and(i2<=i3)) then
          if (d3<>d4)or((d3=d4)and(i3<=i4)) then
            begin
              ans := xx(xx(xx(xx(xxxx(A,B,C,D),E),F),G),Z);
              if ans[0,1] <> 0 then
                begin
                  AA[0][0,1] := AA[0][0,1]+1;
                  AA[ AA[0][0,1] ] := ans;
writeln(name,'"((((xxxx(',T[i1],',',T[i2],',',T[i3],',',T[i4],')',T[i5],')',T[i6],')',T[i7],')',T[i8],')",');
                end;
            end;
       end;         
end;

(*---------------------------------------------------------------------*)












procedure eoneabcdx(d1,d2,d3,d4,d5:longint);

var A,B,C,D,E,ans:poly;
    i1,i2,i3,i4,i5:longint;

begin
     A[0,1] := 1; B[0,1] := 1; C[0,1] := 1; D[0,1] := 1; E[0,1] := 1;
     A[1,1] := 1; B[1,1] := 1; C[1,1] := 1; D[1,1] := 1; E[1,1] := 1;
     for i1 := start[d1] to stop[d1] do
     for i2 := start[d2] to stop[d2] do
     for i3 := start[d3] to stop[d3] do
     for i4 := start[d4] to stop[d4] do
     for i5 := start[d5] to stop[d5] do
       begin
          A[1,2] := i1;
          B[1,2] := i2;
          C[1,2] := i3;
          D[1,2] := i4;
          E[1,2] := i5;
          if (d1<>d3)or((d1=d3)and(i1<i3)) then
          if (d3<>d4)or((d3=d4)and(i3<i4)) then
            begin
              ans := xx(eone(A,B,C,D),E);
              if ans[0,1] <> 0 then
                begin
                  AA[0][0,1] := AA[0][0,1]+1;
                  AA[ AA[0][0,1] ] := ans;
writeln(name,'"(eone(',T[i1],',',T[i2],',',T[i3],',',T[i4],')',T[i5],')",');
                end;
            end;
       end;         
end;

procedure eoneabcdxy(d1,d2,d3,d4,d5,d6:longint);

var A,B,C,D,E,F,ans:poly;
    i1,i2,i3,i4,i5,i6:longint;

begin
     A[0,1] := 1; B[0,1] := 1; C[0,1] := 1; D[0,1] := 1; E[0,1] := 1; F[0,1] := 1;
     A[1,1] := 1; B[1,1] := 1; C[1,1] := 1; D[1,1] := 1; E[1,1] := 1; F[1,1] := 1;
     for i1 := start[d1] to stop[d1] do
     for i2 := start[d2] to stop[d2] do
     for i3 := start[d3] to stop[d3] do
     for i4 := start[d4] to stop[d4] do
     for i5 := start[d5] to stop[d5] do
     for i6 := start[d6] to stop[d6] do
       begin
          A[1,2] := i1;
          B[1,2] := i2;
          C[1,2] := i3;
          D[1,2] := i4;
          E[1,2] := i5;
          F[1,2] := i6;
          if (d1<>d3)or((d1=d3)and(i1<i3)) then
          if (d3<>d4)or((d3=d4)and(i3<i4)) then
            begin
              ans := xx(xx(eone(A,B,C,D),E),F);
              if ans[0,1] <> 0 then
                begin
                  AA[0][0,1] := AA[0][0,1]+1;
                  AA[ AA[0][0,1] ] := ans;
writeln(name,'"((eone(',T[i1],',',T[i2],',',T[i3],',',T[i4],')',T[i5],')',T[i6],')",');
                end;
            end;
       end;         
end;


procedure eoneabcdxyz(d1,d2,d3,d4,d5,d6,d7:longint);

var A,B,C,D,E,F,G,ans:poly;
    i1,i2,i3,i4,i5,i6,i7:longint;

begin
     A[0,1]:=1;B[0,1]:=1;C[0,1]:=1;D[0,1]:=1;E[0,1]:=1;F[0,1]:=1;G[0,1]:=1;
     A[1,1]:=1;B[1,1]:=1;C[1,1]:=1;D[1,1]:=1;E[1,1]:=1;F[1,1]:=1;G[1,1]:=1;
     for i1 := start[d1] to stop[d1] do
     for i2 := start[d2] to stop[d2] do
     for i3 := start[d3] to stop[d3] do
     for i4 := start[d4] to stop[d4] do
     for i5 := start[d5] to stop[d5] do
     for i6 := start[d6] to stop[d6] do
     for i7 := start[d7] to stop[d7] do
       begin
          A[1,2] := i1;
          B[1,2] := i2;
          C[1,2] := i3;
          D[1,2] := i4;
          E[1,2] := i5;
          F[1,2] := i6;
          G[1,2] := i7;
          if (d1<>d3)or((d1=d3)and(i1<i3)) then
          if (d3<>d4)or((d3=d4)and(i3<i4)) then
            begin
              ans := xx(xx(xx(eone(A,B,C,D),E),F),G);
              if ans[0,1] <> 0 then
                begin
                  AA[0][0,1] := AA[0][0,1]+1;
                  AA[ AA[0][0,1] ] := ans;
writeln(name,'"(((eone(',T[i1],',',T[i2],',',T[i3],',',T[i4],')',T[i5],')',T[i6],')',T[i7],')",');
                end;
            end;
       end;         
end;


procedure eoneabcdxyzw(d1,d2,d3,d4,d5,d6,d7,d8:longint);

var A,B,C,D,E,F,G,Z,ans:poly;
    i1,i2,i3,i4,i5,i6,i7,i8:longint;

begin
     A[0,1]:=1;B[0,1]:=1;C[0,1]:=1;D[0,1]:=1;E[0,1]:=1;F[0,1]:=1;G[0,1]:=1;Z[0,1]:=1;
     A[1,1]:=1;B[1,1]:=1;C[1,1]:=1;D[1,1]:=1;E[1,1]:=1;F[1,1]:=1;G[1,1]:=1;Z[1,1]:=1;
     for i1 := start[d1] to stop[d1] do
     for i2 := start[d2] to stop[d2] do
     for i3 := start[d3] to stop[d3] do
     for i4 := start[d4] to stop[d4] do
     for i5 := start[d5] to stop[d5] do
     for i6 := start[d6] to stop[d6] do
     for i7 := start[d7] to stop[d7] do
     for i8 := start[d8] to stop[d8] do
       begin
          A[1,2] := i1;
          B[1,2] := i2;
          C[1,2] := i3;
          D[1,2] := i4;
          E[1,2] := i5;
          F[1,2] := i6;
          G[1,2] := i7;
          Z[1,2] := i7;
          if (d1<>d3)or((d1=d3)and(i1<i3)) then
          if (d3<>d4)or((d3=d4)and(i3<i4)) then
            begin
              ans := xx(xx(xx(xx(eone(A,B,C,D),E),F),G),Z);
              if ans[0,1] <> 0 then
                begin
                  AA[0][0,1] := AA[0][0,1]+1;
                  AA[ AA[0][0,1] ] := ans;
writeln(name,'"((((eone(',T[i1],',',T[i2],',',T[i3],',',T[i4],')',T[i5],')',T[i6],')',T[i7],')',T[i8],')",');

                end;
            end;
       end;         
end;

function ctoe(nn:longint):poly;

var ans:poly;

begin
  if combasis[nn,2] = 0 then
     begin
       ans[0,1] := 1;
       ans[1,1] := 1;
       ans[1,2] := nn;
     end;
  if combasis[nn,2] <> 0 then
     begin
       ans := xx(ctoe(combasis[nn,2]),ctoe(combasis[nn,3]));
     end;
  ctoe := ans;  
end;


function eonetocommutative(n:longint):longint;

var nout,n1,n2,n3,i:longint;

begin
       if basis[n,2] = 0 then
          begin
            nout := n;
          end;
       if basis[n,2] <> 0 then
          begin
            n2 := eonetocommutative(basis[n,2]);
            n3 := eonetocommutative(basis[n,3]);
            nout := 0; 
            for n1 := 1 to combasislength do
                if (combasis[n1,2] = n2) and (combasis[n1,3] = n3) then
                     nout := n1;
           end;
     if nout = 0 then
        begin
          writeln(' in eonetocommutative we have a problen with ');
          writeln(' n is ',n,'   ',T[n]);
          readln(i);
          writeln(' got to use this somehow ',i);
        end; 
     eonetocommutative := nout;
end;

function commutative(n2,n3:longint):longint;

var n1,ans:longint;

begin
   ans := 0;
   for n1:=1 to combasislength do
     if ((combasis[n1,2] = n2) and (combasis[n1,3] = n3))or
        ((combasis[n1,2] = n3) and (combasis[n1,3] = n2)) then
        ans := n1;
   if ans = 0 then
       begin
         writeln(' n2 and n3 ', n2,' ',n3, ' have no commutative product ');
         writeln(' aw shit ');
         writeln(' comT n2 ',comT[n2]);
         writeln(' comT n3 ',comT[n3]);
         readln(ans);
       end;
   commutative := ans;
end;

(* a,b,c,d are commutative terms                                       *)
(* creates ((ab)c)d in commutative terms                               *)
(* looks up the reduction in commutative terms                         *)
(* which should be reduced                                             *)


function ooo(a,b,c,d:longint):poly;
var n,i:longint;
    ans:poly;

begin
     n := commutative(commutative(commutative(a,b),c),d); 
     ans := ctoe(etermfile[n]);
     for i:=1 to ans[0,1] do
       begin
         ans[i,1] := (ans[i,1]*coefficient[n]) mod PRIME;  
       end;
     ooo := ans;
end;

function ooox(a,b,c,d,e:longint):poly;
var n,i:longint;
    ans:poly;

begin
     n := commutative(commutative(commutative(commutative(a,b),c),d),e); 
     ans := ctoe(etermfile[n]);
     for i:=1 to ans[0,1] do
       begin
        ans[i,1] := (ans[i,1]*coefficient[n]) mod PRIME;  
       end;
     ooox := ans;
end;

function oooxy(a,b,c,d,e,f:longint):poly;
var n,i:longint;
    ans:poly;

begin
     n := commutative(commutative(commutative(commutative(commutative(a,b),c),d),e),f); 
     ans := ctoe(etermfile[n]);
     for i:=1 to ans[0,1] do
       begin
        ans[i,1] := (ans[i,1]*coefficient[n]) mod PRIME;  
       end;
     oooxy := ans;
end;

function oooxyz(a,b,c,d,e,f,g:longint):poly;
var n,i:longint;
    ans:poly;

begin
     n := commutative(commutative(commutative(commutative(commutative(commutative(a,b),c),d),e),f),g); 
     ans := ctoe(etermfile[n]);
     for i:=1 to ans[0,1] do
       begin
        ans[i,1] := (ans[i,1]*coefficient[n]) mod PRIME;  
       end;
     oooxyz := ans;
end;

function oooxyzw(a,b,c,d,e,f,g,h:longint):poly;
var n,i:longint;
    ans:poly;

begin
n := commutative(commutative(commutative(commutative(commutative(commutative(commutative(a,b),c),d),e),f),g),h); 
     ans := ctoe(etermfile[n]);
     for i:=1 to ans[0,1] do
       begin
        ans[i,1] := (ans[i,1]*coefficient[n]) mod PRIME;  
       end;
     oooxyzw := ans;
end;


function oooxyzwu(a,b,c,d,e,f,g,h,k:longint):poly;
var n,i:longint;
    ans:poly;

begin
     n := commutative(a,b); 
     n := commutative(n,c); 
     n := commutative(n,d); 
     n := commutative(n,e); 
     n := commutative(n,f); 
     n := commutative(n,g); 
     n := commutative(n,h); 
     n := commutative(n,k); 

     ans := ctoe(etermfile[n]);
     for i:=1 to ans[0,1] do
       begin
         ans[i,1] := (ans[i,1]*coefficient[n]) mod PRIME;  
       end;
     oooxyzwu := ans;
end;


procedure oooabcd(d1,d2,d3,d4:longint);

var ans:poly;
    i1,i2,i3,i4:longint;
    X1,X2,X3,X4,X5,X6:poly;

begin
     for i1 := comstart[d1] to comstop[d1] do
     for i2 := comstart[d2] to comstop[d2] do
     for i3 := comstart[d3] to comstop[d3] do
     for i4 := comstart[d4] to comstop[d4] do
       begin
         if comparts[i1,1]+comparts[i2,1]+comparts[i3,1]+comparts[i4,1] = 1 then
         if comparts[i1,2]+comparts[i2,2]+comparts[i3,2]+comparts[i4,2] = 1 then
         if comparts[i1,3]+comparts[i2,3]+comparts[i3,3]+comparts[i4,3] = 1 then
         if comparts[i1,4]+comparts[i2,4]+comparts[i3,4]+comparts[i4,4] = 1 then
         if comparts[i1,5]+comparts[i2,5]+comparts[i3,5]+comparts[i4,5] <= 5 then
           begin
             X1 := ooo(i1,i2,i3,i4); 
             X2 := ooo(i1,i3,i4,i2); 
             X3 := ooo(i1,i4,i2,i3); 
             X4 := ooo(i1,i3,i2,i4); 
             X5 := ooo(i1,i2,i4,i3); 
             X6 := ooo(i1,i4,i3,i2); 
             ans := X1;
             ans := xApyB(1,ans,+1,X2);
             ans := xApyB(1,ans,+1,X3);
             ans := xApyB(1,ans,-1,X4);
             ans := xApyB(1,ans,-1,X5);
             ans := xApyB(1,ans,-1,X6);
             if ans[0,1] <> 0 then
                  begin
                    AA[0][0,1] := AA[0][0,1]+1;
                    AA[ AA[0][0,1] ] := ans;
 writeln(name,'"EONE(',comT[i1],',',comT[i2],',',comT[i3],',',comT[i4],')",');
                  end;
             end;
      end;
end;

procedure oooabcdx(d1,d2,d3,d4,d5:longint);

var ans:poly;
    i1,i2,i3,i4,i5:longint;
    X1,X2,X3,X4,X5,X6:poly;

begin
     for i1 := comstart[d1] to comstop[d1] do
     for i2 := comstart[d2] to comstop[d2] do
     for i3 := comstart[d3] to comstop[d3] do
     for i4 := comstart[d4] to comstop[d4] do
     for i5 := comstart[d5] to comstop[d5] do
       begin
 if comparts[i1,1]+comparts[i2,1]+comparts[i3,1]+comparts[i4,1]+comparts[i5,1] = 1 then
 if comparts[i1,2]+comparts[i2,2]+comparts[i3,2]+comparts[i4,2]+comparts[i5,2] = 1 then
 if comparts[i1,3]+comparts[i2,3]+comparts[i3,3]+comparts[i4,3]+comparts[i5,3] = 1 then
 if comparts[i1,4]+comparts[i2,4]+comparts[i3,4]+comparts[i4,4]+comparts[i5,4] = 1 then
 if comparts[i1,5]+comparts[i2,5]+comparts[i3,5]+comparts[i4,5]+comparts[i5,5] <= 5 then
           begin
             X1 := ooox(i1,i2,i3,i4,i5); 
             X2 := ooox(i1,i3,i4,i2,i5); 
             X3 := ooox(i1,i4,i2,i3,i5); 
             X4 := ooox(i1,i3,i2,i4,i5); 
             X5 := ooox(i1,i2,i4,i3,i5); 
             X6 := ooox(i1,i4,i3,i2,i5); 
             ans := X1;
             ans := xApyB(1,ans,+1,X2);
             ans := xApyB(1,ans,+1,X3);
             ans := xApyB(1,ans,-1,X4);
             ans := xApyB(1,ans,-1,X5);
             ans := xApyB(1,ans,-1,X6);
             if ans[0,1] <> 0 then
                  begin
                    AA[0][0,1] := AA[0][0,1]+1;
                    AA[ AA[0][0,1] ] := ans;
 writeln(name,'"EONE(',comT[i1],',',comT[i2],',',comT[i3],',',comT[i4],')',comT[i5],'",');
                  end;
             end;
      end;
end;


procedure oooabcdxy(d1,d2,d3,d4,d5,d6:longint);

var ans:poly;
    i1,i2,i3,i4,i5,i6:longint;
    X1,X2,X3,X4,X5,X6:poly;
    i,s1,s2,s3,s4,s5:longint;

begin
     for i1 := comstart[d1] to comstop[d1] do
     for i2 := comstart[d2] to comstop[d2] do
     for i3 := comstart[d3] to comstop[d3] do
     for i4 := comstart[d4] to comstop[d4] do
     for i5 := comstart[d5] to comstop[d5] do
     for i6 := comstart[d6] to comstop[d6] do
       begin
s1 := comparts[i1,1]+comparts[i2,1]+comparts[i3,1]+comparts[i4,1]+comparts[i5,1]+comparts[i6,1];   
s2 := comparts[i1,2]+comparts[i2,2]+comparts[i3,2]+comparts[i4,2]+comparts[i5,2]+comparts[i6,2];   
s3 := comparts[i1,3]+comparts[i2,3]+comparts[i3,3]+comparts[i4,3]+comparts[i5,3]+comparts[i6,3];   
s4 := comparts[i1,4]+comparts[i2,4]+comparts[i3,4]+comparts[i4,4]+comparts[i5,4]+comparts[i6,4];   
s5 := comparts[i1,5]+comparts[i2,5]+comparts[i3,5]+comparts[i4,5]+comparts[i5,5]+comparts[i6,5];   

         if (s1 = 1)and(s2 = 1)and(s3 = 1)and(s4 = 1)and(s5 <= 5) then
           begin
             X1 := oooxy(i1,i2,i3,i4,i5,i6); 
             X2 := oooxy(i1,i3,i4,i2,i5,i6); 
             X3 := oooxy(i1,i4,i2,i3,i5,i6); 
             X4 := oooxy(i1,i3,i2,i4,i5,i6); 
             X5 := oooxy(i1,i2,i4,i3,i5,i6); 
             X6 := oooxy(i1,i4,i3,i2,i5,i6); 
             ans := X1;
             ans := xApyB(1,ans,+1,X2);
             ans := xApyB(1,ans,+1,X3);
             ans := xApyB(1,ans,-1,X4);
             ans := xApyB(1,ans,-1,X5);
             ans := xApyB(1,ans,-1,X6);
             if ans[0,1] <> 0 then
                  begin
                    AA[0][0,1] := AA[0][0,1]+1;
                    AA[ AA[0][0,1] ] := ans;
 writeln(name,'"(EONE(',comT[i1],',',comT[i2],',',comT[i3],',',comT[i4],')',comT[i5],')',comT[i6],'",');
                  end;
             end;
      end;
end;

procedure oooabcdxyz(d1,d2,d3,d4,d5,d6,d7:longint);

var ans:poly;
    i1,i2,i3,i4,i5,i6,i7:longint;
    X1,X2,X3,X4,X5,X6:poly;
    i,s1,s2,s3,s4,s5:longint;

begin
     for i1 := comstart[d1] to comstop[d1] do
     for i2 := comstart[d2] to comstop[d2] do
     for i3 := comstart[d3] to comstop[d3] do
     for i4 := comstart[d4] to comstop[d4] do
     for i5 := comstart[d5] to comstop[d5] do
     for i6 := comstart[d6] to comstop[d6] do
     for i7 := comstart[d7] to comstop[d7] do
       begin
s1 := comparts[i1,1]+comparts[i2,1]+comparts[i3,1]+comparts[i4,1]+comparts[i5,1]+comparts[i6,1]+comparts[i7,1]; 
s2 := comparts[i1,2]+comparts[i2,2]+comparts[i3,2]+comparts[i4,2]+comparts[i5,2]+comparts[i6,2]+comparts[i7,2]; 
s3 := comparts[i1,3]+comparts[i2,3]+comparts[i3,3]+comparts[i4,3]+comparts[i5,3]+comparts[i6,3]+comparts[i7,3]; 
s4 := comparts[i1,4]+comparts[i2,4]+comparts[i3,4]+comparts[i4,4]+comparts[i5,4]+comparts[i6,4]+comparts[i7,4]; 
s5 := comparts[i1,5]+comparts[i2,5]+comparts[i3,5]+comparts[i4,5]+comparts[i5,5]+comparts[i6,5]+comparts[i7,5]; 

         if (s1 = 1)and(s2 = 1)and(s3 = 1)and(s4 = 1)and(s5 <= 5) then
           begin
             X1 := oooxyz(i1,i2,i3,i4,i5,i6,i7); 
             X2 := oooxyz(i1,i3,i4,i2,i5,i6,i7); 
             X3 := oooxyz(i1,i4,i2,i3,i5,i6,i7); 
             X4 := oooxyz(i1,i3,i2,i4,i5,i6,i7); 
             X5 := oooxyz(i1,i2,i4,i3,i5,i6,i7); 
             X6 := oooxyz(i1,i4,i3,i2,i5,i6,i7); 
             ans := X1;
             ans := xApyB(1,ans,+1,X2);
             ans := xApyB(1,ans,+1,X3);
             ans := xApyB(1,ans,-1,X4);
             ans := xApyB(1,ans,-1,X5);
             ans := xApyB(1,ans,-1,X6);
             if ans[0,1] <> 0 then
                  begin
                    AA[0][0,1] := AA[0][0,1]+1;
                    AA[ AA[0][0,1] ] := ans;
 write(name,'"((EONE(',comT[i1],',',comT[i2],',',comT[i3],',',comT[i4],')');
 write(name,comT[i5],')',comT[i6],')',comT[i7],'",');
 writeln(name);
                  end;
             end;
      end;
end;

procedure oooabcdxyzw(d1,d2,d3,d4,d5,d6,d7,d8:longint);

var ans:poly;
    i1,i2,i3,i4,i5,i6,i7,i8:longint;
    X1,X2,X3,X4,X5,X6:poly;
    i,s1,s2,s3,s4,s5:longint;

begin
     for i1 := comstart[d1] to comstop[d1] do
     for i2 := comstart[d2] to comstop[d2] do
     for i3 := comstart[d3] to comstop[d3] do
     for i4 := comstart[d4] to comstop[d4] do
     for i5 := comstart[d5] to comstop[d5] do
     for i6 := comstart[d6] to comstop[d6] do
     for i7 := comstart[d7] to comstop[d7] do
     for i8 := comstart[d8] to comstop[d8] do
       begin
s1 := comparts[i1,1]+comparts[i2,1]+comparts[i3,1]+comparts[i4,1]+comparts[i5,1]+comparts[i6,1]+comparts[i7,1]; 
s1 := s1+comparts[i8,1];

s2 := comparts[i1,2]+comparts[i2,2]+comparts[i3,2]+comparts[i4,2]+comparts[i5,2]+comparts[i6,2]+comparts[i7,2]; 
s2 := s2+comparts[i8,2];

s3 := comparts[i1,3]+comparts[i2,3]+comparts[i3,3]+comparts[i4,3]+comparts[i5,3]+comparts[i6,3]+comparts[i7,3]; 
s3 := s3+comparts[i8,3];

s4 := comparts[i1,4]+comparts[i2,4]+comparts[i3,4]+comparts[i4,4]+comparts[i5,4]+comparts[i6,4]+comparts[i7,4]; 
s4 := s4+comparts[i8,4];

s5 := comparts[i1,5]+comparts[i2,5]+comparts[i3,5]+comparts[i4,5]+comparts[i5,5]+comparts[i6,5]+comparts[i7,5]; 
s5 := s5+comparts[i8,5];

         if (s1 = 1)and(s2 = 1)and(s3 = 1)and(s4 = 1)and(s5 <= 5) then
           begin
             X1 := oooxyzw(i1,i2,i3,i4,i5,i6,i7,i8); 
             X2 := oooxyzw(i1,i3,i4,i2,i5,i6,i7,i8); 
             X3 := oooxyzw(i1,i4,i2,i3,i5,i6,i7,i8); 
             X4 := oooxyzw(i1,i3,i2,i4,i5,i6,i7,i8); 
             X5 := oooxyzw(i1,i2,i4,i3,i5,i6,i7,i8); 
             X6 := oooxyzw(i1,i4,i3,i2,i5,i6,i7,i8); 
             ans := X1;
             ans := xApyB(1,ans,+1,X2);
             ans := xApyB(1,ans,+1,X3);
             ans := xApyB(1,ans,-1,X4);
             ans := xApyB(1,ans,-1,X5);
             ans := xApyB(1,ans,-1,X6);
             if ans[0,1] <> 0 then
                  begin
                    AA[0][0,1] := AA[0][0,1]+1;
                    AA[ AA[0][0,1] ] := ans;
 write(name,'"(((EONE(',comT[i1],',',comT[i2],',',comT[i3],',',comT[i4],')');
 write(name,comT[i5],')',comT[i6],')',comT[i7],')',comT[i8],'",');
 writeln(name);
                  end;
             end;
      end;
end;


procedure oooabcdxyzwu(d1,d2,d3,d4,d5,d6,d7,d8,d9:longint);

var ans:poly;
    i1,i2,i3,i4,i5,i6,i7,i8,i9:longint;
    X1,X2,X3,X4,X5,X6:poly;
    i,s1,s2,s3,s4,s5:longint;

begin
     for i1 := comstart[d1] to comstop[d1] do
     for i2 := comstart[d2] to comstop[d2] do
     for i3 := comstart[d3] to comstop[d3] do
     for i4 := comstart[d4] to comstop[d4] do
     for i5 := comstart[d5] to comstop[d5] do
     for i6 := comstart[d6] to comstop[d6] do
     for i7 := comstart[d7] to comstop[d7] do
     for i8 := comstart[d8] to comstop[d8] do
     for i9 := comstart[d9] to comstop[d9] do
       begin
s1 := comparts[i1,1]+comparts[i2,1]+comparts[i3,1]+comparts[i4,1]+comparts[i5,1]+comparts[i6,1]+comparts[i7,1]; 
s1 := s1+comparts[i8,1]+comparts[i9,1];

s2 := comparts[i1,2]+comparts[i2,2]+comparts[i3,2]+comparts[i4,2]+comparts[i5,2]+comparts[i6,2]+comparts[i7,2]; 
s2 := s2+comparts[i8,2]+comparts[i9,2];

s3 := comparts[i1,3]+comparts[i2,3]+comparts[i3,3]+comparts[i4,3]+comparts[i5,3]+comparts[i6,3]+comparts[i7,3]; 
s3 := s3+comparts[i8,3]+comparts[i9,3];

s4 := comparts[i1,4]+comparts[i2,4]+comparts[i3,4]+comparts[i4,4]+comparts[i5,4]+comparts[i6,4]+comparts[i7,4]; 
s4 := s4+comparts[i8,4]+comparts[i9,4];;

s5 := comparts[i1,5]+comparts[i2,5]+comparts[i3,5]+comparts[i4,5]+comparts[i5,5]+comparts[i6,5]+comparts[i7,5]; 
s5 := s5+comparts[i8,5]+comparts[i9,5];

         if (s1 = 1)and(s2 = 1)and(s3 = 1)and(s4 = 1)and(s5 <= 5) then
           begin
             X1 := oooxyzwu(i1,i2,i3,i4,i5,i6,i7,i8,i9); 
             X2 := oooxyzwu(i1,i3,i4,i2,i5,i6,i7,i8,i9); 
             X3 := oooxyzwu(i1,i4,i2,i3,i5,i6,i7,i8,i9); 
             X4 := oooxyzwu(i1,i3,i2,i4,i5,i6,i7,i8,i9); 
             X5 := oooxyzwu(i1,i2,i4,i3,i5,i6,i7,i8,i9); 
             X6 := oooxyzwu(i1,i4,i3,i2,i5,i6,i7,i8,i9); 
             ans := X1;
             ans := xApyB(1,ans,+1,X2);
             ans := xApyB(1,ans,+1,X3);
             ans := xApyB(1,ans,-1,X4);
             ans := xApyB(1,ans,-1,X5);
             ans := xApyB(1,ans,-1,X6);
             if ans[0,1] <> 0 then
                  begin
                    AA[0][0,1] := AA[0][0,1]+1;
                    AA[ AA[0][0,1] ] := ans;
 write(name,'"((((EONE(',comT[i1],',',comT[i2],',',comT[i3],',',comT[i4],')');
 write(name,comT[i5],')',comT[i6],')',comT[i7],')',comT[i8],')',comT[i9],'",');
 writeln(name);
                  end;
             end;
      end;
end;

 
(*---------------------------------------------------------------------------*)

begin
   combasislength  :=  55653;
   eonebasislength :=  14502;

   assign(filein,'termfile_commutative');
   reset(filein);
   for i:=1 to combasislength do
     begin
       readln(filein,comT[i]);
     end;
   close(filein);

   assign(filein,'COEFFICIENT');
   reset(filein);
   for i:=1 to combasislength do
      begin
        readln(filein,coefficient[i]);
      end;
   close(filein);

   assign(filein,'ETERMFILE');
   reset(filein);
   for i:=1 to combasislength do
      begin
        readln(filein,etermfile[i]);
      end;
   close(filein);

   writeln(' hells bells ');

   assign(filein,'basis_commutative_abcd5e');
   reset(filein);
   readln(filein);
   for i:=1 to combasislength do
     begin
        readln(filein,xreal,combasis[i,2],combasis[i,3],combasis[i,4]);
        combasis[i,1] := round(xreal);
     end;
   close(filein);

   for i:=1 to combasislength do
     begin
       for j:=1 to 5 do
          comparts[i,j] := 0;

       if combasis[i,2] = 0 then
          comparts[i,i] := 1
       else
         for j:=1 to 5 do
           comparts[i,j] := comparts[combasis[i,2],j]+comparts[combasis[i,3],j];
     end;

   for i:=1 to combasislength do
     begin
       if combasis[i,2] = 0 then
         comlog[i] := 1
       else
         comlog[i] := comlog[combasis[i,2]]+comlog[combasis[i,3]];
     end;

   for i:=1 to combasislength do
     comstop[comlog[i]] := i;

   for i:= combasislength downto 1 do
     comstart[comlog[i]] := i;

for i:=1 to 9 do
   begin
     writeln(' comstart and comstop ',comstart[i],' ',comstop[i]);
   end;


   assign(filein,'basis_eone_abcd5e');
   reset(filein);
   readln(filein);
   for i:=1 to eonebasislength do
     begin
       readln(filein,xreal,basis[i,2],basis[i,3],basis[i,4]);
       basis[i,1] := round(xreal);  
     end;

   for i:=1 to eonebasislength do
     begin
       if basis[i,2] = 0 then
         log[i] := 1
       else
         log[i] := log[basis[i,2]]+log[basis[i,3]];
     end;
   
   for i:=1 to eonebasislength do
        stop[log[i]] := i;

   for i:=eonebasislength downto 1 do
        start[log[i]] := i;

   for i:=1 to 9 do
     writeln(i,' ',start[i],' ',stop[i]);

          
   assign(filein,'termfile');
   reset(filein);
   for i:=1 to  eonebasislength do
     begin
        readln(filein,T[i]);
     end;
   close(filein);

   writeln('ONE YOU GREAT DOG ');

   assign(filein,'STRIP_INDEX');
   reset(filein);
   readln(filein,Index[0,1],Index[0,2],Index[0,3],Index[0,4]);
   for i:=1 to Index[0,1] do
     readln(filein,n,index[i,1],index[i,2],index[i,3],index[i,4]);
   close(filein);

   writeln(' the length of index is ',index[0,1]:10);

   writeln('TWO YOU GREAT DOG');
   assign(filein,'STRIP_TERMS');
   reset(filein);
   readln(filein,H[0,1],H[0,2]);
   for i:=1 to H[0,1] do
      readln(filein,n,H[i,1],H[i,2]);
   close(filein); 
   writeln(' the length of H is ',H[0,1]:10);

   writeln('THREE YOU GREAT DOG');

   assign(name,'name');
   rewrite(name);
   writeln(name,'name = {');

   AA[0][0,1] := 0;

(*  Hello  *)

writeln(' JESUS IS COMING ');


for i1 := 1 to 6 do
for i2 := 1 to 7-i1 do
for i3 := 1 to 8-i1-i2 do
for i4 := 1 to 9-i1-i2-i3 do
  if (i2 >= i3)and(i3 >= i4) then
     begin
       oooabcd(i1,i2,i3,i4)
     end;

for i1 := 1 to 5 do
for i2 := 1 to 6-i1 do
for i3 := 1 to 7-i1-i2 do
for i4 := 1 to 8-i1-i2-i3 do
for i5 := 1 to 9-i1-i2-i3-i4 do
  if (i2 >= i3)and(i3 >= i4) then
     begin
       oooabcdx(i1,i2,i3,i4,i5)
     end;

for i1 := 1 to 4 do
for i2 := 1 to 5-i1 do
for i3 := 1 to 6-i1-i2 do
for i4 := 1 to 7-i1-i2-i3 do
for i5 := 1 to 8-i1-i2-i3-i4 do
for i6 := 1 to 9-i1-i2-i3-i4-i5 do
  if (i2 >= i3)and(i3 >= i4) then
     begin
       oooabcdxy(i1,i2,i3,i4,i5,i6);
     end;

for i1 := 1 to 3 do
for i2 := 1 to 4-i1 do
for i3 := 1 to 5-i1-i2 do
for i4 := 1 to 6-i1-i2-i3 do
for i5 := 1 to 7-i1-i2-i3-i4 do
for i6 := 1 to 8-i1-i2-i3-i4-i5 do
for i7 := 1 to 9-i1-i2-i3-i4-i5-i6 do
  if (i2 >= i3)and(i3 >= i4) then
     begin
       oooabcdxyz(i1,i2,i3,i4,i5,i6,i7);
     end;

for i1 := 1 to 2 do
for i2 := 1 to 3-i1 do
for i3 := 1 to 4-i1-i2 do
for i4 := 1 to 5-i1-i2-i3 do
for i5 := 1 to 6-i1-i2-i3-i4 do
for i6 := 1 to 7-i1-i2-i3-i4-i5 do
for i7 := 1 to 8-i1-i2-i3-i4-i5-i6 do
for i8 := 1 to 9-i1-i2-i3-i4-i5-i6-i7 do
  if (i2 >= i3)and(i3 >= i4) then
     begin
       oooabcdxyzw(i1,i2,i3,i4,i5,i6,i7,i8);
     end;


oooabcdxyzwu(1,1,1,1,1,1,1,1,1);

close(name);

(* WOOF *)

(*
const c1 = 2;
      c2 = 3;
      c3 = 7;
      c4 = 7;

 _                       _  _                                     _
|                         ||                                       |
|(c2-1)(c3-1)(c2-c3) (bc) || -p(2p-1) (ad) -(p-1) (ad)e + ((ad)e)e | = 0
|_                       _||_   -1            0          1        _|

*)
A[0,1] := 1; B[0,1] := 1; C[0,1] := 1; D[0,1] := 1; E[0,1] := 1;
A[1,1] := 1; B[1,1] := 1; C[1,1] := 1; D[1,1] := 1; E[1,1] := 1;
A[1,2] := 1; B[1,2] := 2; C[1,2] := 3; D[1,2] := 4; E[1,2] := 5;
    X1 := xx(A,B);
    X2 := xx(C,D);
    X3 := xx(X2,E);
    X4 := xx(X3,E);

    X5 := xApyB( -c3*(2*c3-1),X2,1,xApyB(-(c3-1), X3, +1,X4));
     
    AA[0][0,1] := AA[0][0,1]+1;
    AA[ AA[0][0,1] ] := xx(X1,X5);

seepoly(AA[ AA[0][0,1] ]);

writeln(' IT IS BETTER THAN ROAST TANG');

col := stop[9]; 
writeln(' the length of col is ',col);


    writeln(' HELLO TANG'); 
  assign(fileout,'BUILD_OUT');
  rewrite(fileout);
  writeln(fileout,'A = {');

     for n:=1 to AA[0][0,1] do
        begin
          for j:=1 to col do 
             mat[j] := 0;
       
          for i:=1 to AA[n][0,1] do
            begin
              j := AA[n][i,2];
              mat[j] := mat[j]+AA[n][i,1];
            end;

           if n < AA[0][0,1] then 
             begin
               writeln(fileout,'{');
               for j:=1 to col-1 do
                 writeln(fileout,mat[j]:3,',');
               writeln(fileout,mat[col]:3,'},');
             end;

           if n = AA[0][0,1] then 
             begin
               writeln(fileout,'{');
               for j:=1 to col-1 do
                 writeln(fileout,mat[j]:3,',');
               writeln(fileout,mat[col]:3,'}};');
             end;
        end; 

     writeln(' GOOD TO GO  and col is ',col);

      close(fileout);

end. 
     


