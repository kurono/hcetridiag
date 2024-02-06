/* Main application to solve for 1D transient heat conduction in a rod.

The PDE is discretized using the finite difference algorithm.

The resulting system of linear equations is solved exactly 
using the double sweep tridiagonal matrix method.

The temperature profile along the rod at different time steps 
is plotted as a series of lines. 
The color of each line corresponds to the simulation time, ranging from blue (time = 0) to red (time = max).

(c) Tsivilskiy Ilya, 10.08.2018 */

var RainbowColormap = (function() {
    var _sharp = 10;
    var _power = 2;
    var _ff = 255;

    function RainbowColormap() {

    }

    function gaussian(value, offset) {
       return parseInt(_ff / (1 + Math.pow((value - offset), _power) * _sharp));
    }

    function red(value) {
        return gaussian(value, 1);
    }
    
    function green(value) {
        return gaussian(value, .5);
    }
    
    function blue(value) {
        return gaussian(value, 0);
    }
    
    RainbowColormap.prototype.rgb = function(value) {
        var c = {r: parseInt(red(value)), g: parseInt(green(value)), b: parseInt(blue(value))};
        var cDecimal = ( c.r << 16 ) | ( c.g << 8 ) | c.b;
        var sHex = cDecimal.toString(16);
        //var cHex = parseInt(sHex, 16);
        return sHex;
    };

    return RainbowColormap;
})();

var Dispatcher = (function(){
    var _events = {};

    function Dispatcher() {

    }

    Dispatcher.prototype.addEvent = function(type, callback) {
        var event = document.createEvent('Event');
        event.initEvent(type, true, true);

        document.addEventListener(type, function() {
            callback();
        }, false);

        _events[type] = event;
    };

    Dispatcher.prototype.callEvent = function(type) {
        var event = _events[type];
        document.dispatchEvent(event);
    };

    return Dispatcher;
})();

var Vector2 = (function(){
    var _that;

    function Vector2(x, y) {
        _that = this;

        _that.x = x;
        _that.y = y;
    }

    return Vector2;
})();

var Renderer2d = (function() {
    var _that;
    var _ctx;
    var WIDTH;
    var HEIGHT;

    function init(canvasId) {
        var canvas = document.getElementById(canvasId);

        WIDTH = parseFloat(canvas.width);
        HEIGHT = parseFloat(canvas.height);

        _ctx = canvas.getContext("2d");
    }

    function Renderer2d(canvasId) {
        _that = this;
        init(canvasId);
    }

    Object.defineProperty(Renderer2d.prototype, "width", {
        get: function() {
            return WIDTH;
        }
    });

    Object.defineProperty(Renderer2d.prototype, "height", {
        get: function() {
            return HEIGHT;
        }
    });

    Renderer2d.prototype.clear = function() {
        _ctx.beginPath();
        _ctx.rect(0, 0, WIDTH, HEIGHT);
        _ctx.fillStyle = "#ffffff";
        _ctx.fill();
        _ctx.closePath();
        _ctx.clearRect(0, 0, WIDTH, HEIGHT);
    };

    /*Renderer2d.prototype.screenToCartesian = function(vec2) {
        return new Vector2(vec2.x - WIDTH/2, HEIGHT/2  - vec2.y);
    };

    Renderer2d.prototype.cartesianToScreen = function(vec2) {
        return new Vector2(vec2.x + WIDTH / 2, HEIGHT / 2 - vec2.y);
    };*/

    Renderer2d.prototype.setLineStyle = function(width, color) {
        _ctx.lineWidth = width;
        _ctx.strokeStyle = color;
    };

    Renderer2d.prototype.moveTo = function(vec2) {
        _ctx.moveTo(vec2.x, vec2.y);
    };

    Renderer2d.prototype.lineTo = function(vec2) {
        _ctx.lineTo(vec2.x, vec2.y);
    };

    Renderer2d.prototype.startDrawing = function() {
        _ctx.beginPath();
    };

    Renderer2d.prototype.finishDrawing = function() {
        _ctx.stroke();
    };

    return Renderer2d;
})();

var Main = (function(){
    var _N;
    var _t_start;
    var _L, _k, _rho, _c;
    var _h, _tau;
    var _T = [], _P = [], _Q = [];
    var _time;
    var _that;

    function zeros(n) {
        var a = [];
        var i;
        for (i = 0; i < n; i++) {
            a.push(0);
        }
        return a;
    }

    function initEventDispatcher() {
        _that.dispatcher = new Dispatcher();
        /*_that.dispatcher.addEvent("event1", function() {alert(1)});
         _that.dispatcher.callEvent("event1");*/
    }

    function Main() {
        _that = this;

        initEventDispatcher();
    }

    Object.defineProperty(Main.prototype, "N", {
        get: function() {
            return _N;
        }
    });

    Object.defineProperty(Main.prototype, "T", {
        get: function() {
            return _T;
        }
    });

    Object.defineProperty(Main.prototype, "percentDone", {
        get: function() {
            return _time/_that.t_end;
        }
    });

    Main.prototype.setDefaults = function() {
        _L       = .1;       // rod length [m]
        _k  	 = 410;      // thermal conductivity (Cu)
        _rho     = 8920;     // density
        _c       = 385;      // heat capacity

        _that.T0      = 300;      // initial T [K]
        _that.Tl      = 400;      // T at the left boundary, x=0
        _that.Tr      = 600;      // T at the right boundary, x=_L

        _t_start = 0;		 	  // initial time
        _that.t_end   = 30;       // stop time

        _N       = 30;      // spatial resolution of an 1D mesh
    };

    Main.prototype.initData = function() {
        _h       = _L /(_N-1);    // mesh element size
        _tau     = (_that.t_end - _t_start) / (_N-1);  // time step

        _T = zeros(_N);
        _P = zeros(_N);
        _Q = zeros(_N);
    };

    Main.prototype.solve = function(log) {
        var i;
        var ai, bi, ci, di; // coefficients of canonical form of the tridiagonal SOLE
		var lambda = _k/(_c*_rho)*_tau/(_h*_h);

        // initial T
        for (i = 0; i < _N-1; i++) {
            _T[i] = _that.T0;
        }

        // numerical integration of a transient heat conduction equation by the double-sweep method
        _time = _t_start;

        while (_time < _that.t_end) {
            _time += _tau;
            // compute initial values of sweep coefficients from the left BC
            _P[0] = 0;
            _Q[0] = _that.Tl;

            // compute the rest of sweep coefficients
            for (i = 1; i < _N; i++) {
                // coefficients of canonical form of the tridiagonal SOLE
                ai = lambda;
                bi = - (1 + 2*lambda);
                ci = lambda;
                di = -1 * _T[i];
                // compute sweep coefficients _P[i] and _Q[i]
                _P[i] = - ci / (ai * _P[i-1] + bi);
                _Q[i] = (di - ai * _Q[i-1]) / (ai * _P[i-1] + bi);
            }
            // set the right-side T
            _T[_N] = _that.Tr;
            // compute unknown T by back substitution
            for (i = _N - 1;  i > 0;  i--) {
                _T[i] = _P[i] * _T[i+1] + _Q[i];
            }

            // send global event to make sure the solution has been completed at the current time step
            _that.dispatcher.callEvent("TIMESTEP_DONE");
        }

        if (log === true) {
            alert("t=" + _time.toFixed(1) + " s, T0=" + _T[0].toFixed(1) + " K, T= " +  _T[_N-1].toFixed(1) + " K");
        }
    };

    return Main;
})();

/// Main script ///

var renderer;
var app;

function setAppParameterFromInput(id) {
    app[id] = parseFloat(document.getElementById(id).value);
}

function init() {
    renderer = new Renderer2d("canvas");
    renderer.clear();
    renderer.setLineStyle(2, "#00ccff");

    app = new Main();
    app.setDefaults();
    app.initData();

    var cmap = new RainbowColormap();

    var plotCurrentSolution = function() {
        //console.log("drawPoint");

        // compute T limits
        var minT = Math.min(app.Tl, app.T0, app.Tr);
        var maxT = Math.max(app.Tl, app.T0, app.Tr);

        // plot T
        renderer.startDrawing();
        renderer.setLineStyle(2, "#" + cmap.rgb(app.percentDone));

        var i;
        var x, y;
        var xy;
        for (i = 1; i < app.N; i++) {
            x = i * renderer.width / app.N;
            y = (1 - (app.T[i] - minT) / (maxT - minT)) * renderer.height;
            xy = new Vector2(x, y);

            if (i == 1) renderer.moveTo(xy);
            renderer.lineTo(xy);
        }
        renderer.finishDrawing();
    };

    app.dispatcher.addEvent("TIMESTEP_DONE", plotCurrentSolution);

    app.solve(false);
}

function onSimulateClick() {
    renderer.clear();
    app.setDefaults();
    app.initData();

    setAppParameterFromInput("T0");
    setAppParameterFromInput("Tl");
    setAppParameterFromInput("Tr");

    console.log(app.Tl + " " + app.T0 + " " + app.Tr);

    app.solve(false);
}