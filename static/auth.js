let client = new jso.JSO({
    providerID: "HBP",
    client_id: "af59e7dc-9ed7-4e7e-baac-9ce5cc0fed78",
    redirect_uri: "https://bspg.pa.ibf.cnr.it:14999/model/",
    authorization: "https://services.humanbrainproject.eu/oidc/authorize",
});

function getUserInfo(pagename){
    // reauthorize if needed
    authorize();
    var current_token = sessionStorage.getItem('token');

    fetch('https://services.humanbrainproject.eu/idm/v1/api/user/me',{
            method: 'get',
            headers: {
                "Authorization": "Bearer " + sessionStorage.getItem('token')
            }
    })
    .then(response => response.json())
        .then(data => {
            fetch('/model/log_user/' + data["id"] + '/' + pagename)
        });
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
