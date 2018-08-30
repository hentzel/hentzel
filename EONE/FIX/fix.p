program fix(input,output);

const c1 = 2;
      c2 = 3;
      c3 = 7;
      c4 = 7;

var i,j,n,nn:longint;
    i1,i2,i3,i4:longint;
    xreal:double;
    start,stop:array[1..17] of longint;
    log:array[1..60000] of longint;
    basis:array[1..60000,1..4] of longint;
    filein,fileout,coefficient,etermfile:text;
    coef:array[1..60000,1..16] of longint;
    T:array[1..60000] of string[60];

procedure look(n:longint);

var i1,i2,i3:longint;

begin
   if basis[n,1] < 6 then
     begin
       write(chr(96+basis[n,1]));
     end;

   if basis[n,1] >= 6 then 
      begin
       i2 := basis[n,2];
       i3 := basis[n,3];
       If i2 < i3 then
         begin
           writeln(' shit in display');
           writeln(i2,' ',i3);
           writeln(' paused ');
           readln(i1);
           writeln(' I had to use i1 somehow ',i1);
         end;

       if (i2<6)and(i3<6) then
         write('(',chr(96+i2),chr(96+i3),')');

       if (i2<6)and(i3>=6) then
          begin
            write('(',chr(96+i2));
            look(i3);
            write(')');
          end;

       if (i2>=6)and(i3<6) then
          begin
             write('(');
             look(i2);
             write(chr(96+i3));
             write(')');
          end;

       if (i2>=6)and(i3>=6) then
          begin
            write('(');
            look(i2);
            look(i3);
            write(')');
          end;
     end;
end;
        

begin
   assign(filein,'termfile_commutative');
   reset(filein);
   for i:=1 to  55653  do
     readln(filein,T[i]);
   close(filein);

   assign(filein,'basis_commutative_abcd5e');
   reset(filein);
   readln(filein);
   for i:=1 to   55653 do
     begin
       readln(filein,xreal,i2,i3,i4);
       i1 := round(xreal);
       basis[i,1] := i1;
       basis[i,2] := i2;
       basis[i,3] := i3;
       basis[i,4] := i4;

       If i2 = 0 then log[i] := 1;
       if i2 <> 0 then log[i] := log[i2]+log[i3];
       stop[log[i]] := i;
     end;

   start[1] := 1;
   for i:=1 to 8  do
     start[i+1] := stop[i]+1;

   for i:=1 to   55653  do
      for j:=1 to 16 do
         coef[i,j] := 0; 

(* Hello *)
writeln(' (((((ba)e)c)d)e) ');

   for i:=6 to stop[9] do
     begin
        i2 := basis[basis[i,2],1];
        i3 := basis[basis[i,3],1];

        If (i2<i3) then
          begin
            n  := i2; 
            i2 := i3;
            i3 :=  n;
          end;

        if (i2 = 5)and(i3 = 5) then
          begin
            basis[i,1] := 5;
            basis[i,2] := 0;
            basis[i,3] := 0;
          end;

        for n := 1 to 4 do
           if (i2 = 5)and(i3 = n) then
              begin
                for j:=1 to 16 do 
                   coef[i,j] := coef[basis[i,2],j]+coef[basis[i,3],j];
                coef[i,n] := coef[i,n]+1;
                basis[i,1] := n; 
                basis[i,2] := 0;
                basis[i,3] := 0;
              end;

        if (i2 = 16)and(i3 = 5) then
          begin
            for j:=1 to 16 do 
               coef[i,j] := coef[basis[i,2],j]+coef[basis[i,3],j];
            coef[i,5] := coef[i,5]+1;
            basis[i,1] := i2; 
            basis[i,2] := basis[basis[i2,2],1];
            basis[i,3] := basis[basis[i2,3],1];
          end;

        if (i2 > 6)and(i3 = 5)and(i2<>16) then
           begin
             for j:=1 to 16 do 
               coef[i,j] := coef[basis[i,2],j]+coef[basis[i,3],j];
             basis[i,1] := i;
             for n:=i downto 1 do
               begin
    if ((basis[n,2] = i2)and(basis[n,3] = i3))or((basis[n,2] = i3)and(basis[n,3] = i2)) then
                   basis[i,1] := n;
                 if basis[i2,1] >= basis[i3,1] then
                   begin
                     basis[i,2] := basis[i2,1];
                     basis[i,3] := basis[i3,1];
                   end
                 else
                   begin
                     basis[i,2] := basis[i3,1];
                     basis[i,3] := basis[i2,1];
                   end
               end;
           end;
          

        if (i2 <> 5)and(i3 <> 5)and(i>=7) then
          begin
            for j:=1 to 16 do
               coef[i,j] := coef[basis[i,2],j]+coef[basis[i,3],j];
            basis[i,1] := i; 
            for n:=i downto 1 do
              if ((basis[n,2] = i2)and(basis[n,3] = i3))or
                 ((basis[n,2] = i3)and(basis[n,3] = i2)) then
                basis[i,1] := n;
            if basis[i2,1] >= basis[i3,1] then
              begin
                 basis[i,2] := basis[i2,1];
                 basis[i,3] := basis[i3,1];
              end 
            else
              begin
                basis[i,2] := basis[i3,1];
                basis[i,3] := basis[i2,1];
              end 
          end;
     end;

     assign(coefficient,'COEFFICIENT');
     rewrite(coefficient);

     assign(etermfile,'ETERMFILE');
     rewrite(etermfile);

     for i:=1 to stop[9] do
       begin
         write(i:10,'  ',T[i],'           ');
         write(T[basis[i,1]]);
         for n:=1 to 5 do
           for j:=1 to coef[i,n] do
             write(' c',n:1); 
         write('    ',coef[i,1],' ',coef[i,2],' ',coef[i,3],' ',coef[i,4],' ',coef[i,5]);
         writeln;
       end;

     for i:=1 to stop[9] do
       begin
         writeln(etermfile,basis[i,1]);
         nn := 1;
         for n:=1 to coef[i,1] do
            nn := c1*nn;

         for n:=1 to coef[i,2] do
            nn := c2*nn;

         for n:=1 to coef[i,3] do
            nn := c3*nn;

         for n:=1 to coef[i,4] do
            nn := c4*nn;

         for n:=1 to coef[i,5] do
            nn := (c1+c2-1)*nn;


         writeln(coefficient,nn:1); 

       end;

  close(coefficient);
  close(etermfile);
end.

