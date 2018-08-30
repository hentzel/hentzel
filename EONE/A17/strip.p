program strip(input,output);

type pair = array[1..2] of longint;

     pairlist = array[0..1000] of pair;

     line = array[1..80] of char;

var    

    w : line; 
    i,j:longint;  
    filein,fileout:text;
    product,dependence:boolean;
    H:array[0..250000000] of pair;
    index:array[0..3000000,1..4] of longint;    (* b1*b2 start stop *)

function atoi(x:char):longint;
  
   begin
      atoi := 0;
      if x = '0' then atoi := 0;
      if x = '1' then atoi := 1;
      if x = '2' then atoi := 2;
      if x = '3' then atoi := 3;
      if x = '4' then atoi := 4;
      if x = '5' then atoi := 5;
      if x = '6' then atoi := 6;
      if x = '7' then atoi := 7;
      if x = '8' then atoi := 8;
      if x = '9' then atoi := 9;
   end;

procedure star(w:line);

var i,n,a,b:longint;
    pq:array[1..4] of longint; 
  
begin 
    n := 0; 
    for i:=1 to 60 do
       begin
         if (w[i] = '(') or (w[i] = ')') then
            begin
               n := n+1;
               pq[n] := i;
            end;
       end; 
    a := 0;
    for i:=pq[1]+2 to pq[2]-1 do
       begin
         a := 10*a + atoi(w[i]); 
       end;
    b := 0;
    for i:=pq[3]+2 to pq[4]-1 do
       begin
         b := 10*b + atoi(w[i]); 
       end;
   index[ index[0,1], 4 ] := H[0,1]; 
   index[0,1] := index[0,1]+1;
   index[ index[0,1], 1 ] := a;
   index[ index[0,1], 2 ] := b;
   index[ index[0,1], 3 ] := H[0,1]+1;
   index[ index[0,1], 4 ] := H[0,1]+1;  (* <=== Will be overwritten *)
end;

procedure tail(w:line);

var n:longint;
    a:longint;   
    slot:longint;

begin
    slot := 1;
    a := 0;
    for n:=1 to 60 do
       begin
         if (w[n] = 'b') or (w[n] = '+') then
            begin
              if slot = 1 then
                begin
                  H[0,1] := H[0,1]+1; 
                  H[ H[0,1],1 ] := a;    
                  slot := 2;
                end
              else
                  begin
                    H[ H[0,1],2] := a; 
                    slot := 1;
                  end;
              a := 0;
            end;
         if (w[n] <> '+' ) and (w[n] <> ' ') and (w[n] <> 'b')  then
              a := 10*a + atoi(w[n]);
       end;
     If slot = 2 then
        H[ H[0,1],2 ] := a;
end;  



begin
    H[0,1] := 0;
    index[0,1] := 0;
    assign(filein,'mult_eone_17a');
    reset(filein);
    for i:=1 to  4315273   do
       begin
 writeln(i);
         product := false;
         dependence := true;

         for j:=1 to 60 do
            w[j] := ' ';

         j := 0;
         while not eoln(filein) do
            begin
              j:=j+1;
              read(filein,w[j]);
            end;
          readln(filein);
         
         for j:=1 to 60 do
           begin
             if w[j] = '*' then
                begin
                   dependence := false;
                   product := true;
                end;
           end;
             if dependence then
                 tail(w);
             if product then
                   star(w);
       
     end;
     close(filein);

writeln('TANG IS GOOD');

   assign(fileout,'STRIP_TERMS');
   rewrite(fileout);
   writeln(fileout,H[0,1]:10,H[0,2]:10);
   for i:=1 to H[0,1] do
      writeln(fileout,i:3,' ',H[i,1]:5,' ',H[i,2]:5);
   close(fileout);

   assign(fileout,'STRIP_INDEX');
   rewrite(fileout);
   writeln(fileout,index[0,1]:10,index[0,2]:10,index[0,3]:10,index[0,4]:10);
   for i:=1 to index[0,1] do
      writeln(fileout,i:3,' ',index[i,1],' ',index[i,2],' ',index[i,3],' ',index[i,4]);
   close(fileout); 

end.    

