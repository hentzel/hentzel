program build(input,output);

const PRIME = 251; (* NOT USED *)

type poly = array[0..5000,1..2] of longint;       (* coefficient basis     *)
                                                   (*  p[0,1] is the length *)
var
    i,j,k,n:longint;    
    H:array[0..140000000,1..2] of longint;
             (*  139799736   *)
    index:array[0.. 503721 ,1..4] of longint;
    name,filein,fileout:text;
    A,B,C,D,E,F,G,Z,ans:poly;     
    AA:array[0..1000] of poly;

    mat:array[1..1000,1..10000] of longint;  (* #columns is 7293*)
    col:longint;
    i1,i2,i3,i4,i5,i6,i7,i8:longint;
    x1,x2,x3,x4,x5,x6,x7,x8:char;
    w:array[1..8] of longint;
    start,stop:array[1..17] of longint;
    T:array[1..251889] of string[60];
    xreal:double;
    degone,degtwo:longint;

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

function simplify(C:poly):poly;

var i,j:longint;
    D:poly;

begin
    for i:=1 to C[0,1]-1 do
      for j:=i+1 to C[0,1] do
        if C[i,2] = C[j,2] then
           begin
             C[i,1] := C[i,1]+C[j,1];
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

var n,i,z:longint;
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
                ans[ ans[0,1] ,1] := C[k,1]*A[i,1]*B[j,1]; 
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
procedure eoneabcd(d1,d2,d3,d4:longint);

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
          if (d1<>d3)or((d1=d3)and(i1<i3)) then
          if (d3<>d4)or((d3=d4)and(i3<i4)) then
            begin
              ans := eone(A,B,C,D);
              if ans[0,1] <> 0 then
                begin
                  AA[0][0,1] := AA[0][0,1]+1;
                  AA[ AA[0][0,1] ] := ans;
writeln(name,'"eone(',T[i1],',',T[i2],',',T[i3],',',T[i4],')",');
                end;
            end;
       end;         
end;


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

 
(*---------------------------------------------------------------------------*)

begin
   assign(filein,'basis_eone_17a');
   reset(filein);
   readln(filein);
   for i:=1 to  6418 do
     begin
       readln(filein,xreal,i2,i3,i4);
       stop[i4] := i;  
     end;
   start[1] := 1;
   for i:=1 to 17 do
     start[i+1] := stop[i]+1;

   for i:=1 to 17 do
     writeln(i,' ',start[i],' ',stop[i]);

   writeln(' Paused for santa claus');
   readln(i);

   assign(filein,'termfile');
   reset(filein);
   for i:=1 to  stop[17] do
     readln(filein,T[i]);
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

writeln('THREE YOU GREAT DOG');
assign(name,'name');
rewrite(name);
writeln(name,'name = {');
AA[0][0,1] := 0;
writeln(' the length of H is ',H[0,1]:10);

(*  Hello  *)

        A[0,1] := 1;
        A[1,1] := 1;
        A[1,2] := 1;

        B[0,1] := 1;
        B[1,1] := 1;
        B[1,2] := 2;

        C[0,1] := 1;
        C[1,1] := 1;
        C[1,2] := 3;

        D[0,1] := 1;
        D[1,1] := 1;
        D[1,2] := 4;

        E[0,1] := 1;
        E[1,1] := 1;
        E[1,2] := 5;
 
        F[0,1] := 1;
        F[1,1] := 1;
        F[1,2] := 6;

        G[0,1] := 1;
        G[1,1] := 1;
        G[1,2] := 7;

        Z[0,1] := 1; 
        Z[1,1] := 1;
        Z[1,2] := 8;

(* puke *)
(* 2x12 is zero *)
(* 3x11 is zero *)
(* 4x10 is zero *)
(* 5x9  is zero *)
(* 6x8  is one  *)

(* 2x13         *)
(* 3x12         *)
(* 4x11 is zero        *)
(* 5x10 is zero        *)
(* 6x9  three of them  *)
(* 7x8          *)

(* 4x12 is zero *)
(* 5x10 is zero *)
(* 5x11 is one  but does not factor *)
(* 5x12 is ????????????????         *)
(* 6x7 is zero *)
(* 6x8 is one  *)
(* 6x9 is three but just one six-factor  *)
(* 6x10  ????????????                    *)
(* 7x8 is one  *)
(* 7x9 is three *)

degone := 7; 
degtwo := 8;

AA[0][0,1] := 0;
   begin
         for i1 := start[degone] to stop[degone] do
         for i2 := start[degtwo] to stop[degtwo] do
           begin 
             F[1,2] := i1;
             G[1,2] := i2;
             Z := xx(F,G); 
writeln(' santa ',Z[0,1]);
             AA[0][0,1] := AA[0][0,1] + 1;
             AA[ AA[0][0,1] ] := Z;
             If(i1+i2 < stop[degone]+stop[degtwo]) then
                writeln(name,'"',T[i1],'"*"',T[i2],'",')
             else
                writeln(name,'"',T[i1],'"*"',T[i2],'"};')
           end;
   end;     


close(name);


For i:=1 to AA[0][0,1] do
   begin
     writeln(i:3,' ',AA[i][0,1]);
   end;

Writeln(' paused showing the lengths of the expansions.. BLOW ');
readln(i);

writeln(' IT IS BETTER THAN ROAST TANG');

col := stop[degone+degtwo]-stop[degone+degtwo-1];
writeln(' the length of col is ',col);

     for i:=1 to AA[0][0,1] do
         for j:=1 to col do 
           mat[i,j] := 0;

    writeln(' HELLO TANG'); 

     for n:=1 to AA[0][0,1] do
        begin
          for i:=1 to AA[n][0,1] do
            begin
              j := AA[n][i,2]-stop[degone+degtwo-1];

  writeln('fuck i and j ',i,' ',j);

              mat[n,j] := mat[n,j]+AA[n][i,1];
    
            end;
        end; 

     writeln(' GOOD TO GO  and col is ',col);


(*
              for i:=1 to col do
         begin
       for n:=1 to AA[0][0,1] do
              write(mat[n,i]:3);
           writeln;
         end;       
*)

     assign(fileout,'BUILD_OUT');
     rewrite(fileout);
     writeln(fileout,'A = {');
     for i:=1 to AA[0][0,1]-1 do
         begin
           writeln(fileout,'{');
           for j:=1 to col-1 do
              writeln(fileout,mat[i,j]:3,',');
           writeln(fileout,mat[i,col]:3,'},');
         end; 
       
      writeln(fileout,'{');
      for j:=1 to col-1 do
        writeln(fileout,mat[AA[0][0,1],j]:3,',');
      writeln(fileout,mat[AA[0][0,1],col]:3,'}};');
      close(fileout);

for i:=1 to AA[0][0,1] do
  begin
     displaypoly(AA[i]);
    writeln('==============================================');
  end;
end. 
     


