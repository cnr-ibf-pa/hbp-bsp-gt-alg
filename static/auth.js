var client_id = "";
var redirect_uri = "";

if (url.includes("https://bsp-gtalg.cineca.it")){
	console.log("Loading setting client id PROD (vm on CINECA)");
    client_id = '9788129e-d0a1-4e1a-bb17-6eb07b779609';
    redirect_uri= "https://bsp-gtalg.cineca.it/model/";
} else if (url.includes("https://bspg.pa.ibf.cnr.it:14999")){
	console.log("Loading client id DEV (vm on bspg.pa.ibf.cnr.it)");
    client_id = "af59e7dc-9ed7-4e7e-baac-9ce5cc0fed78";
    redirect_uri = "https://bspg.pa.ibf.cnr.it:14999/model/";
} 

let client = new jso.JSO({
    providerID: "HBP",
    client_id : client_id,
    redirect_uri : redirect_uri,
    authorization: "https://services.humanbrainproject.eu/oidc/authorize",
});

function getUserInfo(pagename){
    // reauthorize if needed
    authorize();
    var current_token = sessionStorage.getItem('token');
    var userid = "";

    fetch('https://services.humanbrainproject.eu/idm/v1/api/user/me',{
            method: 'get',
            headers: {
                "Authorization": "Bearer " + sessionStorage.getItem('token')
            }
    })
    .then(response => response.json())
        .then(data => {
            fetch('/model/log_user/' + data["id"] + '/' + pagename)
            userid = data["id"]
        });
    return userid
}

function authorize(){
    try {
        client.callback();
    } catch (e) {
        console.warn("Issue decoding the token");
    }

    var authorization = client.getToken();
    authorization.then((session) => {
        var access_token = session.access_token;
        sessionStorage.setItem('token', access_token)
    });
}
