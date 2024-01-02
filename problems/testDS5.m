% show distribution of feasbile and infeasible solutions
prob = DS5(5, 4);

xu = linspace(1, 2, 100);
xl1 = linspace(0, 1, 100);
xlrest = zeros(100, 9);

xl = [xl1', xlrest];
xu = xu';

XU = repelem(xu, 100, 1);
XL = repmat(xl, 100, 1);

[f, c] =  prob.evaluate_u(XU, XL);

feasible = c <= 0;
infeasible = c>0;

ff = f(feasible, :);
fi = f(infeasible, :);

scatter(ff(:, 1), ff(:, 2), 30, "red", 'filled');
hold on;
scatter(fi(:, 1), fi(:, 2), 30, 'green', 'filled');
