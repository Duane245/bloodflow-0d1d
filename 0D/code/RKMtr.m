function y = RKMtr(F_xy,x,y,h,A,b)
k_1 = F_xy(x,y,A,b);
 
k_2 = F_xy(x+0.5*h,y+0.5*h*k_1,A,b);
 
k_3 = F_xy((x+0.5*h),(y+0.5*h*k_2),A,b);
 
k_4 = F_xy((x+h),(y+k_3*h),A,b);
 
 
y = y + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h; 
end
