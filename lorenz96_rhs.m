function dydt = lorenz96_rhs(~,y)
	dydt = (circshift(y, -1) - circshift(y, 2)).*circshift(y,1) - y +8;
end
