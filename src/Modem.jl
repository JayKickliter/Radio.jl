################################################################################
#                           _  _ ____ ___  ____ _  _                           #
#                           |\/| |  | |  \ |___ |\/|                           #
#                           |  | |__| |__/ |___ |  |                           #
################################################################################                                     

type Modem
    modulation::Modulation
    samplesPerSymbol::Real
    txFilter::FIRFilter
    rxFilter::FIRFilter
end


function Modem( modulation::Modulation; samplesPerSymbol::Real = 2, pulseShape = rrcos, excessBandwidth = 0.35 )
    samplesPerSymbol >= 2 || error( "samplesPerSymbol must be greater than 2" )
    resampRate = isinteger( samplesPerSymbol ) ? int( samplesPerSymbol )//1 : samplesPerSymbol
    numTaps    = 32 * 10 * samplesPerSymbol
    span       = int( numTaps/samplesPerSymbol )
    taps       = pulseShape( excessBandwidth, span, samplesPerSymbol )
    txFilter   = FIRFilter( taps, resampRate )  
    rxFilter   = FIRFilter( taps, 1/resampRate )  
    
    Modem( modulation, samplesPerSymbol, txFilter, rxFilter )
end

function modulate( modem::Modem, data::AbstractVector )
    symbols       = modulate( modem.modulation, data )
    resampSymbols = filt( modem.txFilter, symbols )
end

function demodulate( modem::Modem, symbols::AbstractVector )
    resampSymbols = filt( modem.rxFilter, symbols )
    data          = demodulate( modem.modulation, resampSymbols )
end
