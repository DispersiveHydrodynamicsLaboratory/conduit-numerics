% Make a smooth box
function[f] = make_smooth_box(vstretch,ampl,boxmax,x0)

f = @(x) ampl/2 * tanh((x-x0+boxmax/2)/vstretch) - ampl/2*tanh((x-x0-boxmax/2)/vstretch);
