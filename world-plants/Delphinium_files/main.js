

var lastOrder = false;
$(function() {
    $('.lazy').Lazy();
});

function initLazyLoad() {

    $('.lazy').Lazy();

}
function toggleShow(id) {
    $('#'+id).toggle();
}

// Smooth scroll for the menu and links with .scrollto classes
function activeHandler() {
$('.nav-menu a, #mobile-nav a, .scrollto').on('click', function() {

    if (location.pathname.replace(/^\//, '') == this.pathname.replace(/^\//, '') && location.hostname == this.hostname) {
        var target = $(this.hash);
        console.log(target);
        //$(target).addClass('highlight');
        if (target.length) {
            var top_space = 330;

            if ($('#header').length) {
                top_space = $('#header').outerHeight();

                if (! $('#header').hasClass('header-scrolled')) {
                    top_space = top_space - 20;
                }
            }

            $('html, body').animate({
                scrollTop: target.offset().top - top_space
            });

            if ($(this).parents('.nav-menu').length) {
                $('.nav-menu .menu-active').removeClass('menu-active');
                $(this).closest('li').addClass('menu-active');
            }

            if ($('body').hasClass('mobile-nav-active')) {
                $('body').removeClass('mobile-nav-active');
                $('#mobile-nav-toggle i').toggleClass('fa-times fa-bars');
                $('#mobile-body-overly').fadeOut();
            }
            return false;
        }
    }
});
}

var listLength = 0;
$(function() {

    //$("#selectBox2").val(3);
    console.log('On');
    console.log($("#selectBox2 option").length);
    listLength = $("#selectBox2 option").length;
});
function doScrollTop() {
    $('html, body').animate({scrollTop : 0},20, 'easeInOutExpo');
    return false;
}


function initSelectBox() {
    listLength = $("#selectBox2 option").length;
}


function checkButtons() {
    var lastOption;
    var firstOption;
    let curVal = $("#selectBox2").val();
    let i = 0;
    $("#selectBox2 > option").each(function() {
        i++;
        if(i == 2) firstOption = this.value;
        lastOption = this.value;
    });

    if(curVal == lastOption) {
        console.log('is Last');
        $("#nextBtn").addClass('isDisabled');
        $("#backBtn").removeClass('isDisabled');
    } else if(curVal == firstOption || curVal == '') {
        console.log('is first');
        $("#backBtn").addClass('isDisabled');
        $("#nextBtn").removeClass('isDisabled');
    } else {
        $("#nextBtn, #backBtn").removeClass('isDisabled');
    }
    console.log('lastOption:'+lastOption+' cur'+curVal);

}


function nextOptionV2() {
    console.log('------------------');
    clearSearch();

    $(".backBtn").removeClass('isDisabled');
    //$("#select").val('eins');
    //$("#select").change();
    //console.log('df'+ $("#selectBox2").val());
    var endReached = false;
    var orderLength = $("#order option").length;
    var familyLength = $("#family option").length;
    var genusLength = $("#genus option").length;
    var iG = 0;
    var iF = 0;
    var stop = false;
    //var trigger;
    var triggerCurGenus;
    var triggerCurFamily;
    var familyValue = $("#family").val();
    var fldValue = $("#genus").val();

    if(fldValue == null) {
        fldValue = 0;
    }




    $("#genus > option").each(function() {
            iG++;
            //console.log(genusLength + '=='+iG);
            if(triggerCurGenus == true) {
                $("#genus").val(this.value);
                changeFilter('changeGenus');
                return false;
                //console.log(genusLength + '=='+iG);
                //return false;
            }
        //console.log('lastOrder'+lastOrder);
            if (genusLength == iG ) {
                var genusLast = true;
                $("#family > option").each(function () {
                        iF++;

                    //console.log('lastOrder2'+lastOrder);

                    if (familyLength == iF && lastOrder === true) {
                        console.log('END of Family reached IN last Order');
                        //return false;

                    }
                    if (familyLength == (iF+1) && lastOrder === true) {
                      //  console.log('Next is END');
                      //  return false;
                    }

                    if (familyLength == iF) {
                        var familyLast = true;
                        console.log('fire NEXTGENUS->Order');
                        if(endReached === false) {
                            changeFilter('nextGenus');
                        }
                    }

                    if(triggerCurFamily == true) {
                        //$("#family").val(this.value);
                        //$("#family").trigger("chosen:updated");
                        console.log('fire NEXTGENUS');
                        if(endReached === false) {
                            changeFilter('nextGenus');
                        }

                        return false;
                    } else {

                        if(familyValue == this.value && familyLength == (iF-1)) {
                            console.log('Next is END');
                        }

                        if(familyValue == this.value) {
                            triggerCurFamily = true;
                            console.log('familyValue '+familyValue);
                        }
                    }

                }
            );
            } else {
                let iF2 = 0;
                $("#family > option").each(function () {
                    iF2++;
                    //console.log('CHECKKKK familyValue:'+familyValue +' --'+this.value);
                    if(familyValue == this.value && familyLength == (iF2) && genusLength == (iG+1) && lastOrder === true) {
                        console.log('THIS IS ABSOLUT END');
                        console.log('NEXT IIIIS END');
                        $(".nextBtn").addClass('isDisabled');
                        endReached = true;
                        return false;
                    } else {
                        endReached = false;
                        $(".nextBtn").removeClass('isDisabled');
                    }
                }
                );




                if(fldValue == this.value) {
                    triggerCurGenus = true;
                    console.log('check '+fldValue);
                }


            }

        }

    );
}

function nextOptionV2_V2() {
    clearSearch();

    $(".backBtn").removeClass('isDisabled');
    //$("#select").val('eins');
    //$("#select").change();
    //console.log('df'+ $("#selectBox2").val());
    var orderLength = $("#order option").length;
    var familyLength = $("#family option").length;
    var genusLength = $("#genus option").length;
    var iG = 0;
    var iF = 0;
    //var trigger;
    var triggerCurGenus;
    var triggerCurFamily;
    var familyValue = $("#family").val();
    var fldValue = $("#genus").val();

    if(fldValue == null) {
        fldValue = 0;
    }




    $("#genus > option").each(function() {
        iG++;
        if (genusLength == iG) {
            var genusLast = true;
            $("#family > option").each(function () {
                    iF++;
                    if (familyLength == iF) {
                        var familyLast = true;
                    }
                    if(triggerCurFamily == true) {
                        $("#family").val(this.value);
                        $("#family").trigger("chosen:updated");
                        changeFilter('changeGenus');
                    }
                    triggerCurFamily = true;
                }
            );
        }

    }

    );
}


function nextOptionV2_OLD() {
    clearSearch();

    $(".backBtn").removeClass('isDisabled');
    //$("#select").val('eins');
    //$("#select").change();
    //console.log('df'+ $("#selectBox2").val());
    var listLength = $("#genus option").length;
    var i = 0;
    var trigger;
    var fldValue = $("#genus").val();
    if(fldValue == null) {
        fldValue = 0;
    }




    $("#genus > option").each(function() {
        i++;
        if(trigger) {
            console.log(listLength +'=='+ i);
            if(listLength == i) {
                //last entry in list
                //$(".nextBtn").addClass('isDisabled');
                console.log('NEXT Off');
                //jaxon_callJaxon('v2changeFilter','nextGenus',jxn.getFormValues('filterForm'));
            }
            console.log('select '+ i);
            $("#genus").val(this.value);
            $("#genus").trigger("chosen:updated");
            //jaxon_callJaxon('v2selectGenus',jxn.getFormValues('filterForm'));

            changeFilter('changeGenus');

            return false;
        }



        if(fldValue == this.value) {
            console.log(this.text + ' ' + this.value + ' << '+$("#genus").val() +'=='+ this.value);
            trigger = true;
        } else {
            console.log(this.text + ' ' + this.value + '  '+$("#genus").val() +'=='+ this.value);
        }
    });
}

function nextOption() {
    clearSearch();
    //$("#backBtn").prop('disabled', false);
    $("#backBtn").removeClass('isDisabled');
    //$("#select").val('eins');
    //$("#select").change();
    var listLength = $("#selectBox2 option").length;


    var i = 0;
    var trigger;

    $("#selectBox2 > option").each(function() {
        i++;
        if(trigger) {
            console.log(listLength +'=='+ i);
            if(listLength == i) {
                //$("#nextBtn").prop('disabled', true);
                $("#nextBtn").addClass('isDisabled');
                console.log('NEXT Off');
            }
            //console.log('select '+ i);
            $("#selectBox2").val(this.value);
            $("#selectBox2").trigger("chosen:updated");

            let form = jxn.getFormValues('newSourceForm');

            let tstamp = Math.floor(Date.now() / 1000);
            let curUrl = window.location.href;
            let parts = curUrl.split('#');

            //let newUrl = parts[0] + '#'+tstamp;
            let newUrl = parts[0] + '#g-'+ form.selectBox2;
            let state = {
                cmd : 'selectFam',
                action : 'normal',
                tabStatus: 'filter-tab',
                //filter: filter,
                form : form
            }
            console.log('state');
            console.log(state);
            history.pushState(state, tstamp, newUrl);
            //jaxon_callJaxon('v2changeFilter', action, state.form);

            jaxon_callJaxon('selectFam',form);
            return false;
        }
        if($("#selectBox2").val() == this.value) {
            console.log(this.text + ' ' + this.value + ' << ');
            trigger = true;
        } else {
            console.log(this.text + ' ' + this.value);
        }
    });
}

function backOption() {
    clearSearch();
    //$("#nextBtn").prop('disabled', false);
    $("#nextBtn").removeClass('isDisabled');

    //$("#select").val(4);

    var before = '';
    var i = 0;

    $("#selectBox2 > option").each(function() {
        i++;
        if($("#selectBox2").val() == this.value) {
            if(i <= 3) {
            //    $("#backBtn").prop('disabled', true);
                $("#backBtn").addClass('isDisabled');
            }
            console.log(this.text + ' ' + this.value + ' << ');
            $("#selectBox2").val(before);
            $("#selectBox2").trigger("chosen:updated");



            let form = jxn.getFormValues('newSourceForm');

            let tstamp = Math.floor(Date.now() / 1000);
            let curUrl = window.location.href;
            let parts = curUrl.split('#');

            //let newUrl = parts[0] + '#'+tstamp;
            let newUrl = parts[0] + '#g-'+ form.selectBox2;
            let state = {
                cmd : 'selectFam',
                action : 'notuse',
                tabStatus: 'filter-tab',
                //filter: filter,
                form : form
            }
            console.log('state');
            console.log(state);
            history.pushState(state, tstamp, newUrl);




            jaxon_callJaxon('selectFam',jxn.getFormValues('newSourceForm'));
            return false;
        } else {
            before = this.value;
            console.log(this.text + ' ' + this.value);
        }
    });
}

function backOptionV2() {
    clearSearch();
    //$("#nextBtn").prop('disabled', false);
    $(".nextBtn").removeClass('isDisabled');

    //$("#select").val(4);

    var before = '';
    var i = 0;


    $("#genus > option").each(function() {
        i++;
        if($("#genus").val() == this.value) {
            if(i <= 3) {
                //    $("#backBtn").prop('disabled', true);
                $(".backBtn").addClass('isDisabled');
            }
            console.log(this.text + ' ' + this.value + ' << ');
            $("#genus").val(before);
            $("#genus").trigger("chosen:updated");
            if(before > 0) {
                //jaxon_callJaxon('v2selectGenus',jxn.getFormValues('filterForm'));
                changeFilter('changeGenus');
            }
            return false;
        } else {
            before = this.value;
            console.log(this.text + ' ' + this.value);
        }
    });
}

function clearSearch() {
    $('#searchFld').val('');
}


lightbox.option({
    'resizeDuration': 100,
    'fadeDuration': 20,
    'imageFadeDuration': 100,
    'wrapAround': true,
    'albumLabel': 'Bild %1 von %2'
})


