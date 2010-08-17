function out = 3DView()
    
    z = circle([0,0], 100, .5);
    profin(:,1:2) = z(:, 1:2);
    profin(:,3) = 0;
    profin(:,4:5) = z(:, 1:2);
    profin(:,6) = 0;
    profiles = profileinterp(profin,15);
    
    
    
end