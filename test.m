fig = uifigure;                              
sld = uislider(fig,'Value', 2,'ValueChanging', @sliderMoving);
sld.Limits = [-10, 10];

function sliderMoving(~, event)
    x = event.Value;
    plot(x, y);
    disp(x);                           
end