function plot_p_curve(p_dist)
nb_bins = 100;

hold on

hist(p_dist, nb_bins);

H = hist(p_dist, nb_bins);
plot([.05 .05], [0 max(H(:))], 'r', 'linewidth', 2)

text(.7, max(H(:)), sprintf('p(p<.05)=%.04f', mean(p_dist<.05)))

set(gca, 'xtick', 0:.1:1,'xticklabel', 0:.1:1)
xlabel('p-value')
axis([0 1 0 max(H(:))*1.1])
end