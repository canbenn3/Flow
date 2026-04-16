// Scale meters to inches for us silly americans who measure rain in the imperial system
double meter_to_inch(double meter)
{
    return meter * 39.3701;
}

double inch_to_meter(double inch)
{
    return inch * .0254;
}
