var head = document.getElementsByTagName('head')[0];
var ga_script_g = document.createElement('SCRIPT');
var ga_script = document.createElement('SCRIPT');
ga_script_g.async = true;
ga_script.type = 'text/javascript';

var url = window.location.href;

if (url.includes("https://bsp-gtalg.cineca.it")){
	console.log("Loading ga for PROD (vm on CINECA)")
    ga_script_g.src = 'https://www.googletagmanager.com/gtag/js?id=UA-91794319-9';
    ga_script.src = '/static/gtag.js';
} else {
    console.log("Not loading google analytics tags: probably you are not running this from the production site");
} 

head.prepend(ga_script);
head.prepend(ga_script_g);
