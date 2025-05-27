function[x, Iter, Info] = jacobi(n, A, b, Toler, IterMax)
  Info = 0;
  x = zeros(n, 1);
  x = transpose(x);
  v = zeros(n, 1);
  for i=1:n
    x(i) = b(i) / A(i,i);
  endfor
  Iter = 0;
  while (1)
    Iter = Iter + 1;
    for i=1:n
      Soma = 0;
      for j=1:n
        if (i != j)
          Soma = Soma + A(i,j) * x(j);
        endif
      endfor
      v(i) = (b(i) - Soma) / A(i,i);
    endfor
    NormaNum = 0;
    NormaDen = 0;
    for i=1:n
      t = abs(v(i) - x(i));
      if (t > NormaNum)
        NormaNum = t;
      endif
      if (abs(v(i)) > NormaDen)
        NormaDen = abs(v(i));
      endif
      x(i) = v(i);
    endfor
    NormaRel = NormaNum / NormaDen;
    NormaRel
    if (NormaRel <= Toler || Iter >= IterMax)
      break;
    endif
  endwhile
  if (NormaRel <= Toler)
    Info = 0;
  else
    Info = 1;
  endif
endfunction

