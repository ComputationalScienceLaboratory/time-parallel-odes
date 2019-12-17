%my own implementation of a finite differences verification routine for the
%gradient
function maxdiff = findiffverify(x, fun, grad)
    h = 1e-8;
    n = numel(x);
    maxdiff = 0;
    A = eye(n); %grab standard basis vectors from here
    gradeval = grad(x);
    diffs = zeros(n, 1);
    parfor i = 1:n %test in all coordinate directions
        stepf = x + h*A(:, i);
        stepb = x - h*A(:, i);
        fd = (fun(stepf) - fun(stepb))/(2*h);
        reldiff = abs(fd - gradeval(i))/abs(fd);
        %disp(i)
        disp(reldiff)
        diffs(i) = reldiff;
    end
    disp("Max relative difference between finite difference approximation and grad function:");
    disp(max(diffs));
    disp("Mean relative difference:");
    disp(mean(diffs));
end
