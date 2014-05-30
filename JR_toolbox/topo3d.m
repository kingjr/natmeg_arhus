function h = topo3d(pnt,val)
% h = topo3d(pnt,val)
% simplified ft_plot_topo3d by JeanRémi King
h = patch('Vertices', pnt, 'Faces', projecttri(pnt, 'delaunay'), 'CData', val, 'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceLighting', 'none');
end