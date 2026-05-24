function dx = FunLeftHeart1(t,y,A,b)
dx = A*[y(1);y(2);y(3);y(4)]+b;
end
